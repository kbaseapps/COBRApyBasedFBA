class FBAPipeline:
    
    def __init__(self):
        self.is_run_fva = False #fva
        self.is_all_reversible = False #all_reversible
        self.is_pfba = False #minimize_flux
        self.is_single_ko = False #simulate_ko
        self.model = None #COBRApy model object
        self.media_supplement_list = []
        self.feature_ko_list = []
        self.reaction_ko_list = []
        self.custom_bound_list = []
        self.target_reaction = "bio1"
        self.output_id = None #id of the returning FBA object (string)
    
    @staticmethod
    def fromKBaseParams(params):
        pipeline = FBAPipeline()
        ## configure method from params
        ### pipeline.is_run_fva = params[...]
        return pipeline
    
    
    def run(self):
        ##media_supplement_list => [],
        
        #Implement knockouts
        ##feature_ko_list => [],
        ##reaction_ko_list => [],
        #Implement custom bounds
        ##custom_bound_list => [],
        
        #If requested, make all reactions reversible
        
        #Set objective
        
        #Optimize and save this objective for FBA solution
        
        if self.is_pfba:
            #Run PFBA
            pass
        else:
            #optmize
            pass
        
        
        if self.is_run_fva:
            #Run FVA
            ##objective_fraction => 0.1,
            pass
        
        if self.is_single_ko:
            #Simulate all single gene knockouts
            pass
        
        #return result

        pass
        
