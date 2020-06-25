import cobra
import cobrakbase

class FBAPipeline:
    
    def __init__(self):
        self.is_run_fva = False #fva
        self.is_all_reversible = False #all_reversible
        self.is_pfba = False #minimize_flux
        self.is_single_ko = False #simulate_ko
        self.model = None #COBRApy model object
        self.fraction_of_optimum = 1.0
        self.media_supplement_list = []
        self.feature_ko_list = []
        self.reaction_ko_list = []
        self.custom_bound_list = []
        self.target_reaction = 'bio1'
        self.output_id = None #id of the returning FBA object (string)
    
    @staticmethod
    def fromKBaseParams(params):
        pipeline = FBAPipeline()

        # Configure method from params
        self.is_run_fva = params['fva']
        self.is_all_reversible = params['all_reversible']
        self.is_pfba = params['minimize_flux']
        self.is_single_ko = params['simulate_ko']
        self.fraction_of_optimum = params['objective_fraction']
        self.media_supplement_list = params['media_supplement_list']
        self.feature_ko_list = params['feature_ko_list']
        self.reaction_ko_list = params['reaction_ko_list']
        self.custom_bound_list = [] # TODO: doesn't seem to be integrated into UI
        self.target_reaction = params['target_reaction']

        # TODO: add uptake member vars?
        #    'max_c_uptake':          0.,
        #    'max_n_uptake':          0.,
        #    'max_p_uptake':          0.,
        #    'max_s_uptake':          0.,
        #    'max_o_uptake':          0.,
        #    'default_max_uptake':    0.

        self.model = None # TODO: figure this out
    
        return pipeline
    
    def run(self):

        if self.model is None:
            raise ValueError('model must be set before calling run')

        # TODO: what data does this contain? What do we do with it?
        # media_supplement_list => [],
        
        # TODO: only do this if is_single_ko?
        # TODO: what is difference between feature_ko and reaction_ko?
        # Knockouts
        # TODO: this impl probably doesn't work
        # cobra.manipulations.delete_model_genes(self.model, self.feature_ko_list)
        # cobra.manipulations.delete_model_genes(self.model, self.reaction_ko_list)
        
        # Add custom bounds to reactions
        for rct_id, lb, ub in self.custom_bound_list:
            # TODO: should we check if rct_id is in the reactions?
            reaction = self.model.reactions.get_by_id(rct_id)
            reaction.lower_bound = lb
            reaction.upper_bound = ub

        # If requested, make all reactions reversible
        if self.is_all_reversible:
            pass # TODO https://cobrapy.readthedocs.io/en/latest/faq.html?highlight=reversibility#How-do-I-change-the-reversibility-of-a-Reaction?
        
        #Set objective
        #converting the input string which is a reaction ID and setting the flux of this reaction as the objective
        #Optimize and save this objective for FBA solution

        # TODO: can is_pfba and is_run_fva both be True?

        if self.is_pfba:
            from cobra.flux_analysis import pfba
            sol = pfba(self.model,
                       fraction_of_optimum=self.fraction_of_optimum,
                       objective=None,
                       reactions=None)

        if self.is_run_fva:
            from cobra.flux_analysis import flux_variability_analysis as fva
            sol = fva(self.model,
                      reaction_list=None,
                      loopless=False,
                      fraction_of_optimum=self.fraction_of_optimum,
                      pfba_factor=None)

        # TODO: add support for cobra.flux_analysis.geometric_fba ?
        
        if self.is_single_ko:
            #Simulate all single gene knockouts
            pass

        # TODO: does self.output_id get set here?
        #return result (KBaseFBABuilder object)
            
        return
        
