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
        # If true, run FVA algorithm on model.
        self.is_run_fva = False

        # If true, make all reactions in model reversible.
        self.is_all_reversible = False

        # If true, run pFBA to compute FBA solution.
        self.is_pfba = False

        # If true, run gene knockout algorithm on each gene
        # to determine the set of the essential genes.
        self.is_single_ko = False

        # If true, run loopless FBA with CycleFreeFlux
        # algorithm inplace of regular FBA. Will not run
        # when is_pfba is true.
        self.is_loopless_fba = False

        # If true, run loopless FVA with CycleFreeFlux.
        self.is_loopless_fva = False

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
        self.max_uptakes = {atom: 0. for atom in self.UPTAKE_ATOMS}

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

        # TODO: write comment
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
        p.is_run_fva = params['fva']
        p.is_pfba = params['minimize_flux']
        p.is_single_ko = params['simulate_ko']
        p.is_loopless_fba = params['loopless_fba']
        p.is_loopless_fva = params['loopless_fva']
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
        p.custom_bound_list = [] # TODO: doesn't seem to be integrated into UI
        # TODO: model editor example tuple bounds input

        return p

    def configure_media(self, model, media):

        # TODO: should we allow user specified default_max_uptake for complete media?
        #       if so, then check if default uptake != 0
        if media.name == 'Complete':
            self.default_max_uptake = 100.

        if not math.isclose(self.default_max_uptake, 0):
            for ex_flux in model.medium:
                model.reactions.get_by_id(ex_flux).lower_bound = -1 * self.default_max_uptake

        constrs = {atom: 0. for atom in self.UPTAKE_ATOMS}

        for rct_id in model.medium:

            ex_rct = model.reactions.get_by_id(rct_id)
            compound = list(ex_rct.metabolites)[0]
            cmp_atoms = atom_count(compound.formula)

            for atom in self.UPTAKE_ATOMS:
                atom_occurences = cmp_atoms.get(atom)
                if atom_occurences:
                    constrs[atom] += ex_rct.flux_expression * atom_occurences

        # TODO: Have discussion about how to set lb for complete media or if the
        #       compound is in a media formulation
        for atom in self.UPTAKE_ATOMS:
            model.add_cons_vars(
                model.problem.Constraint(constrs[atom],
                                         lb=0,
                                         ub=self.max_uptakes[atom]))

    def run(self, model, media):
        """This function mutates model."""

        # TODO: is there a particular order we must do these steps in?
        #       e.g. should we set all reactions reversible as the last step?

        # Select optimization solver
        model.solver = self.solver
        if self.solver == 'coinor_cbc':
            # Use all available processors
            model.solver.configuration.threads = -1

        # Add custom bounds to reactions
        for rct_id, lb, ub in self.custom_bound_list:
            if rct_id in model.reactions:
                rct = model.reactions.get_by_id(rct_id)
                rct.lower_bound, rct.upper_bound = lb, ub

        # If specified, make all reactions reversible
        if self.is_all_reversible:
            for rct in model.reactions:
                rct.upper_bound = self.MAX_BOUND
                rct.lower_bound = self.MAX_BOUND * -1

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

        # TODO: handle media_supplement_list

        # Set objective
        if self.target_reaction:
            model.objective = self.target_reaction

        # Set direction of the objective function, maximize by default.
        if self.is_minimize_objective:
            model.objective.direction = 'min'

        # Compute FBA solution
        if self.is_pfba:
            from cobra.flux_analysis import pfba
            fba_sol = pfba(model,
                           fraction_of_optimum=self.fraction_of_optimum_pfba)
        elif self.is_loopless_fba:
            # Run CycleFreeFlux algorithm
            fba_sol = cobra.flux_analysis.loopless_solution(model)
        else:
            # Run vanilla FBA
            fba_sol = model.optimize()

        # If specified, compute FVA solution
        fva_sol = None
        if self.is_run_fva:
            from cobra.flux_analysis import flux_variability_analysis as fva
            fva_sol = fva(model,
                          loopless=self.is_loopless_fva,
                          fraction_of_optimum=self.fraction_of_optimum_fva)
        
        # If specified, simulate all single gene knockouts
        essential_genes = set()
        if self.is_single_ko:
            essential_genes = cobra.flux_analysis.variability. \
                              find_essential_genes(model, threshold=1e-11)

        # Convert COBRApy model to kbase format
        fba_builder = KBaseFBABuilder.fromCobra(self.output_id,
                                                model,
                                                fba_sol,
                                                media,
                                                self.workspace)

        b = fba_builder.build()
        # TODO: add fva_sol and essential_genes to this object in cobrakbase
        return b
        

if __name__ == '__main__':

    params = {
        'fbamodel_id':           'Lactococcus_lactis_model',
        'fbamodel_workspace':    'chenry:narrative_1574200759205',
        'media_id':              'RichDefinedAnaerobic.media2',
        'media_workspace':       'chenry:narrative_1574200759205',
        'target_reaction':       '',
        'fba_output_id':         None,
        'workspace':             'chenry:narrative_1574200759205',
        'fva':                   False,
        'minimize_flux':         False,
        'simulate_ko':           False,
        'all_reversible':        False,
        'feature_ko_list':       [], #['L192589', 'L182555'],
        'reaction_ko_list':      [], #['rxn05625_c0', 'rxn00545_c0'],
        'media_supplement_list': [],
        'objective_fraction':    1.,
        'max_c_uptake':          0.,
        'max_n_uptake':          0.,
        'max_p_uptake':          0.,
        'max_s_uptake':          0.,
        'max_o_uptake':          0.,
        'default_max_uptake':    0.
    }
