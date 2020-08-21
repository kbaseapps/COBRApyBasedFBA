import math
import cobra
import cobrakbase
from cobrakbase.modelseed.utils import atom_count
from cobrakbase.core.kbase_fba_builder import KBaseFBABuilder

class FBAPipeline:

    # Bound for making reactions reversible
    MAX_BOUND = 1000

    # Tuple of uptake atoms
    UPTAKE_ATOMS = ('C', 'N', 'P', 'S', 'O')
    
    def __init__(self):
        # If true, make all reactions in model reversible.
        self.is_all_reversible = False

        # If true, run gene knockout algorithm on each gene
        # to determine the set of the essential genes.
        self.is_single_ko = False

        # Sets FBA type. Either 'pFBA', 'Loopless FBA' or 'FBA'
        self.fba_type = 'pFBA'

        # Sets FVA type. Either 'FVA', 'Loopless FVA' or 'Neither'
        self.fva_type = 'FVA'

        # If true, the objective will be minimized.
        # Otherwise the objective is maximized.
        self.is_minimize_objective = False

        # For pFBA, Must be <= 1.0. Requires that the objective value
        # is at least the fraction times maximum objective value.
        self.fraction_of_optimum_pfba = 1.0

        # For FVA, Must be <= 1.0. Requires that the objective value
        # is at least the fraction times maximum objective value.
        self.fraction_of_optimum_fva = 0.1

        # Specify optimization solver, options: [coinor_cbc, glpk].
        self.solver = 'coinor_cbc'

        # Specify objective reaction to optimize.
        self.target_reaction = ''

        # Max uptakes for exchange reactions
        self.max_uptakes = {atom: None for atom in self.UPTAKE_ATOMS}

        # Used with complete media. Default is 0, if complete media
        # change from 0 to 100. Otherwise if user specifies
        # default_max_uptake != 0, every exchange flux has a lower
        # bound of negative default_max_uptake.
        self.default_max_uptake = 0.

        # Knocks out specified genes (which could knockout reactions)
        self.feature_ko_list = []

        # Knockout specified reactions.
        self.reaction_ko_list = []

        # Custom bounds to add to reaction in model
        self.custom_bound_list = []

        # Compounds to add to media
        self.media_supplement_list = []

        # Kbase ID of the returning FBA object.
        self.output_id = ''

        # Kbase workspace name where FBA app is run.
        self.workspace = ''

        # NOTE: boolean variables are passed in as ints
        #       with 0,1 value during runtime.

    @staticmethod
    def fromKBaseParams(params):
        p = FBAPipeline()

        # Kbase environment params
        p.output_id = params['fba_output_id']
        p.workspace = params['fbamodel_workspace']

        # Optimization and algorithmic params
        p.solver = params['solver']
        p.fba_type = params['fba_type']
        p.fva_type = params['fva_type']
        p.is_single_ko = params['simulate_ko']
        p.target_reaction = params['target_reaction']
        p.is_all_reversible = params['all_reversible']
        p.is_minimize_objective = params['minimize_objective']
        p.fraction_of_optimum_fva = params['fraction_of_optimum_fva']
        p.fraction_of_optimum_pfba = params['fraction_of_optimum_pfba']

        # Uptakes
        p.max_uptakes['C'] = params['max_c_uptake']
        p.max_uptakes['N'] = params['max_n_uptake']
        p.max_uptakes['P'] = params['max_p_uptake']
        p.max_uptakes['S'] = params['max_s_uptake']
        p.max_uptakes['O'] = params['max_o_uptake']
        p.default_max_uptake = params['default_max_uptake']

        # Check if list params contain data. If so parse, else use default []
        if params['media_supplement_list']:
            p.media_supplement_list = params['media_supplement_list'].split(',')
        if params['reaction_ko_list']:
            p.reaction_ko_list = params['reaction_ko_list'].split(',')
        if params['feature_ko_list']:
            p.feature_ko_list = params['feature_ko_list'].split(',')
        # custom_bound_list input format [{'custom_reaction_id': ['rct1'], 'custom_lb': float, 'custom_ub': float}]
        if params['custom_bound_list']:
            p.custom_bound_list = [(item['custom_reaction_id'][0], item['custom_lb'], item['custom_ub'])
                                    for item in params['custom_bound_list']]

        return p

    def configure_media(self, model, media):

        if media.name == 'Complete':
            self.default_max_uptake = 100.

        if not math.isclose(self.default_max_uptake, 0):
            for ex_flux in model.medium:
                model.reactions.get_by_id(ex_flux).lower_bound = -1 * self.default_max_uptake

        # Only add constraints when user specifies to, otherwise max_uptakes will be None
        constrs = {atom: 0. for atom in self.UPTAKE_ATOMS if self.max_uptakes[atom] is not None}

        # Return early if user doesn't specify any max uptakes
        if not constrs:
            return

        for rct_id in model.medium:

            ex_rct = model.reactions.get_by_id(rct_id)
            compound = list(ex_rct.metabolites)[0]
            cmp_atoms = atom_count(compound.formula)

            for atom in constrs:
                atom_occurences = cmp_atoms.get(atom)
                if atom_occurences:
                    constrs[atom] += ex_rct.reverse_variable * atom_occurences

        for atom in constrs:
            model.add_cons_vars(
                model.problem.Constraint(constrs[atom],
                                         lb=0,
                                         ub=self.max_uptakes[atom]))

    def run(self, model, media):
        """This function mutates model."""

        # Select optimization solver
        model.solver = self.solver
        if self.solver == 'coinor_cbc':
            # Use all available processors
            model.solver.configuration.threads = -1

        # If specified, make all reactions reversible
        if self.is_all_reversible:
            for rct in model.reactions:
                rct.upper_bound = self.MAX_BOUND
                rct.lower_bound = self.MAX_BOUND * -1

        # Add custom bounds to reactions
        for rct_id, lb, ub in self.custom_bound_list:
            if rct_id in model.reactions:
                rct = model.reactions.get_by_id(rct_id)
                rct.lower_bound, rct.upper_bound = lb, ub

        # Filter out user specified genes to ko that are not in model.
        self.feature_ko_list = list(filter(lambda gene: gene in model.genes,
                                           self.feature_ko_list))
        # Knockout specified genes
        if self.feature_ko_list:
            cobra.manipulation.delete_model_genes(model, self.feature_ko_list)

        for rct_id in self.reaction_ko_list:
            if rct_id in model.reactions:
                rct = model.reactions.get_by_id(rct_id)
                rct.lower_bound, rct.upper_bound = 0, 0

        # Update exchange reaction variable bounds based on user
        # specified max uptakes. Add max uptake contraints.
        self.configure_media(model, media)

        # Set objective
        if self.target_reaction:
            model.objective = self.target_reaction

        # Set direction of the objective function, maximize by default.
        if self.is_minimize_objective:
            model.objective.direction = 'min'

        # Compute FBA solution
        if self.fba_type == 'pFBA':
            from cobra.flux_analysis import pfba
            fba_sol = pfba(model,
                           fraction_of_optimum=self.fraction_of_optimum_pfba)
        elif self.fba_type == 'Loopless FBA':
            # Run CycleFreeFlux algorithm
            fba_sol = cobra.flux_analysis.loopless_solution(model)
        else:
            # Run vanilla FBA
            fba_sol = model.optimize()

        # If specified, compute FVA solution
        fva_sol = None
        if self.fva_type != 'Neither':
            from cobra.flux_analysis import flux_variability_analysis as fva
            fva_sol = fva(model,
                          loopless=self.fva_type == 'Loopless FVA',
                          fraction_of_optimum=self.fraction_of_optimum_fva)
        
        # If specified, simulate all single gene knockouts
        essential_genes = set()
        if self.is_single_ko:
            essential_genes = cobra.flux_analysis.variability. \
                              find_essential_genes(model, threshold=1e-11)

        # Convert COBRApy model to kbase format
        fba_builder = KBaseFBABuilder.from_cobra(self.output_id,
                                                 model,
                                                 fba_sol,
                                                 media,
                                                 self.workspace)

        kbase_fba_obj = fba_builder.with_cobra_fva_solution(fva_sol).build()
        
        # TODO: essential_genes to this object in cobrakbase
        return kbase_fba_obj, fva_sol, fba_sol, essential_genes
