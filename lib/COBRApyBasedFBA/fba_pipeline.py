import math
import cobra
import cobrakbase
from cobrakbase.core.kbase_fba_builder import KBaseFBABuilder

class FBAPipeline:

    MAX_BOUND = 1000
    
    def __init__(self):
        self.is_run_fva = False # run fva
        self.is_all_reversible = False # all_reversible
        self.is_pfba = False # minimize_flux
        self.is_single_ko = False #s imulate_ko
        self.fraction_of_optimum = 1.0
        self.media_supplement_list = []
        self.feature_ko_list = [] # This knocks out genes (which could knockout reactions)
        self.reaction_ko_list = [] # This knocks out features
        self.custom_bound_list = []
        self.target_reaction = ''
        self.output_id = '' # id of the returning FBA object
        self.workspace = ''

    @staticmethod
    def fromKBaseParams(params):
        p = FBAPipeline()

        # Configure method from params
        p.is_run_fva = params['fva']
        p.is_all_reversible = params['all_reversible']
        p.is_pfba = params['minimize_flux']
        p.is_single_ko = params['simulate_ko']
        p.fraction_of_optimum = params['objective_fraction']
        p.media_supplement_list = params['media_supplement_list']
        p.feature_ko_list = params['feature_ko_list']
        p.reaction_ko_list = params['reaction_ko_list']
        p.custom_bound_list = [] # TODO: doesn't seem to be integrated into UI
        p.target_reaction = params['target_reaction']
        p.workspace = params['workspace']

        return p

    def run(self, model, media):
        """This function mutates model."""

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

        # Knockouts
        if self.feature_ko_list:
            cobra.manipulation.delete_model_genes(model, self.feature_ko_list)
        for rct_id in self.reaction_ko_list:
            # TODO: should we check if rct_id is in the reactions?
            #       depends on interface is it plain text? then yes
            rct = model.reactions.get_by_id(rct_id)
            rct.lower_bound, rct.upper_bound = 0, 0

        # Set max uptakes
        # TODO: loop through media and set max uptakes
        # TODO: add uptake member vars? Might have to parse smiles. This mutates the model constraints.
        #       data structure of media: https://narrative.kbase.us/#spec/type/KBaseBiochem.Media
        #    'max_c_uptake':          0.,
        #    'max_n_uptake':          0.,
        #    'max_p_uptake':          0.,
        #    'max_s_uptake':          0.,
        #    'max_o_uptake':          0.,
        #    'default_max_uptake':    0.  # TODO: Ask Chris what this is

        # Set objective
        # converting the input string which is a reaction ID and setting the flux of this reaction as the objective
        # Optimize and save this objective for FBA solution
        # TODO: may have bugs. test this
        if self.target_reaction:
            model.objective = self.target_reaction

        if self.is_pfba:
            from cobra.flux_analysis import pfba
            fba_sol = pfba(model,
                           fraction_of_optimum=self.fraction_of_optimum,
                           objective=None, # TODO: do we need these params?
                           reactions=None)
        else:
            # Run vanilla FBA
            fba_sol = model.optimize()

        fva_sol = None
        if self.is_run_fva:
            from cobra.flux_analysis import flux_variability_analysis as fva
            fva_sol = fva(model,
                          fraction_of_optimum=self.fraction_of_optimum)
        
        # If specified, simulate all single gene knockouts
        essential_genes = set()
        if self.is_single_ko:
            for gene in model.genes:
                cobra.manipulation.delete_model_genes(model, [gene])
                sol = model.optimize()
                if sol.status != 'optimal' or math.isclose(sol.objective_value, 0):
                    essential_genes.add(gene.id)
                cobra.manipulation.undelete_model_genes(model)
                # TODO: add a getter function in KBaseFBABuilder def gene_esential(gene): query set -> bool

        # TODO: temporary return for testing
        return fba_sol, fva_sol, essential_genes

        fba_builder = KBaseFBABuilder.fromCobra(self.output_id,
                                                model,
                                                fba_sol,
                                                media,
                                                self.workspace)
        # TODO: add fva_sol and essential_genes to this object in cobrakbase
            
        return fba_builder.build()
        
