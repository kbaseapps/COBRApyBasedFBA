import math
import cobra
import cobrakbase
import pandas as pd
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

        constrs = {atom: 0. for atom in self.UPTAKE_ATOMS}

        for rct_id in model.medium:

            ex_rct = model.reactions.get_by_id(rct_id)
            compound = list(ex_rct.metabolites)[0]
            cmp_atoms = atom_count(compound.formula)

            for atom in self.UPTAKE_ATOMS:
                atom_occurences = cmp_atoms.get(atom)
                if atom_occurences:
                    constrs[atom] += ex_rct.reverse_variable * atom_occurences

        for atom in self.UPTAKE_ATOMS:
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
        
        # TODO: add fva_sol and essential_genes to this object in cobrakbase
        return fba_builder.build(), fva_sol, fba_sol # TODO: Temporary tuple return. add to fva_builder


# TODO: see above TODOs, fva_sol, fba_sol should be part of fba_object
def build_report(pipeline, model, fva_sol, fba_sol):
    """Build output report and return string of html."""
    import os
    import jinja2

    # Helper funcs for formating
    yes_no_format = lambda x: 'Yes' if x else 'No'
    nan_format = lambda x: x if x else 'NaN'

    fba_type = 'pFBA' if pipeline.is_pfba else 'FBA'

    # Seperate exchange reaction ids
    ex_rcts = fva_sol.loc[fva_sol.index.str[:2] == 'EX'].index
    rcts = fva_sol.loc[fva_sol.index.str[:2] != 'EX'].index

    # Select ATP metabolite
    if 'atp_c' in model.metabolites:
        df = model.metabolites.atp_c.summary().to_frame()
    elif 'cpd00002' in model.metabolites:
        df = model.metabolites.cpd00002.summary().to_frame()
    else:
        # Create empty data frame to display void ATP summary
        df = pd.DataFrame()

    atp_summary = []
    for index, row in zip(df.index, df.itertuples()):
        atp_summary.append([*index, *row[1:]])


    # TODO: add correct vals
    fba_model = 'FBA_MODEL'
    media = 'MEDIA'

    context = {'summary':     [x[1:] for x in model.summary().to_frame().itertuples()],
               'atp_summary': atp_summary,
               'overview':    [{'name': 'Model',                'value': fba_model},
                               {'name': 'Media',                'value': media},
                               {'name': 'Optimization status',  'value': fba_sol.status},
                               {'name': 'Objective',            'value': model.objective},
                               {'name': 'Objective value',      'value': fba_sol.objective_value},
                               {'name': 'Number of reactions',  'value': len(fva_sol)},
                               {'name': 'Number of compounds',  'value': len(model.metabolites)},
                               {'name': 'FBA type',             'value': fba_type},
                               {'name': 'Loopless FBA',         'value': yes_no_format(pipeline.is_loopless_fba)},
                               {'name': 'Loopless FVA',         'value': yes_no_format(pipeline.is_loopless_fva)},
                               {'name': 'FVA fraction of opt.', 'value': pipeline.fraction_of_optimum_fva},
                               {'name': 'Reversible reactions', 'value': yes_no_format(pipeline.is_all_reversible)},
                               {'name': 'Single gene KO',       'value': yes_no_format(pipeline.is_single_ko)},
                               {'name': 'Gene KO',              'value': len(pipeline.feature_ko_list)},
                               {'name': 'Reaction KO',          'value': len(pipeline.reaction_ko_list)},
                               {'name': 'Custom bounds',        'value': len(pipeline.custom_bound_list)},
                               {'name': 'Media supplement',     'value': len(pipeline.media_supplement_list)},
                               {'name': 'Solver',               'value': pipeline.solver.upper()}
                               ],
               'reactions': json.dumps([{'id': rct_id,
                                         'min_flux': fva_sol.minimum[rct_id],
                                         'max_flux': fva_sol.maximum[rct_id],
                                         'equation': rct.reaction,
                                         'name': nan_format(rct.name)}
                                       for rct_id, rct in zip(rcts, map(model.reactions.get_by_id, rcts))]),
               'ex_reactions': json.dumps([{'id': rct_id,
                                            'min_flux': fva_sol.minimum[rct_id],
                                            'max_flux': fva_sol.maximum[rct_id],
                                            'equation': rct.reaction,
                                            'name': nan_format(rct.name)}
                                          for rct_id, rct in zip(ex_rcts, map(model.reactions.get_by_id, ex_rcts))]),
           }

    template_dir = os.path.dirname(os.path.realpath(__file__))
    template_file = 'template.html'

    env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(template_dir),
            autoescape=jinja2.select_autoescape(['html', 'xml']))

    # Return string of html
    return env.get_template(template_file).render(context)
