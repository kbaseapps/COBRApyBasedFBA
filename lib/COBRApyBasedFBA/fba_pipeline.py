import math
import cobra
import cobrakbase
from cobrakbase.core.kbase_fba_builder import KBaseFBABuilder

class FBAPipeline:

    # Bound for making reactions reversible
    MAX_BOUND = 1000
    
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

        # Knocks out specified genes (which could knockout reactions)
        self.feature_ko_list = []

        # Knockout specified reactions.
        self.reaction_ko_list = []

        # TODO: write comment
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

        # Check if list params contain data. If so parse, else use default []
        if params['media_supplement_list']:
            p.media_supplement_list = params['media_supplement_list'].split(',')
        if params['reaction_ko_list']:
            p.reaction_ko_list = params['reaction_ko_list'].split(',')
        if params['feature_ko_list']:
            p.feature_ko_list = params['feature_ko_list'].split(',')
        p.custom_bound_list = [] # TODO: doesn't seem to be integrated into UI

        return p

    def run(self, model, media):
        """This function mutates model."""

        # Select optimization solver
        model.solver = self.solver
        if self.solver == 'coinor_cbc':
            # Use all available processors
            model.solver.configuration.threads = -1

        # Add custom bounds to reactions
        for rct_id, lb, ub in self.custom_bound_list:
            # TODO: should we check if rct_id is in the reactions?
            #       depends on interface is it plain text? then yes
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
            # TODO: should we check if rct_id is in the reactions?
            #       depends on interface is it plain text? then yes
            rct = model.reactions.get_by_id(rct_id)
            rct.lower_bound, rct.upper_bound = 0, 0

        # Set max uptakes
        # TODO: add uptake member vars? Might have to parse smiles. This mutates the model constraints.
        #       data structure of media: https://narrative.kbase.us/#spec/type/KBaseBiochem.Media
        #    'max_c_uptake':          0.,
        #    'max_n_uptake':          0.,
        #    'max_p_uptake':          0.,
        #    'max_s_uptake':          0.,
        #    'max_o_uptake':          0.,
        #    'default_max_uptake':    0.   # Used when user picks complete media. default is 0, media tells which compounds can be consumed
        #                                     if complete media change from 0 to 100. if they have default max uptake != 0, every exchange flux
        #                                     has a lower bound that is negative default_max_uptake.

        # if media.name == 'Complete':
        #     # assume media contains everyhting
        #     # set defult_max_uptake to 100

        # if defult_max_uptake != 0:
        #     for exchange_flux in model: # 'EX' prefix
        #         flux.lb = -1* default_max_uptake

        # for c in model.compounds:
        #     if c in uptakes:
        #         # model.media gives compounds for uptakes
        #         c.formula
        #     TODO: check total_carbon sum(exchangeflux * carbon coefficient).
        #           add constr: total_carbon < max_c_uptake

        # TODO: handle media_supplement_list


        # TODO: Look at compounds and check the formula in cobra model. name of field is formula

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
