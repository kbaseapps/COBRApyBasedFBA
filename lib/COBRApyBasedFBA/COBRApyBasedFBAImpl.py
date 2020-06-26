# -*- coding: utf-8 -*-
#BEGIN_HEADER
from COBRApyBasedFBA.fba_pipeline import FBAPipeline
#END_HEADER


class COBRApyBasedFBA:
    '''
    Module Name:
    COBRApyBasedFBA

    Module Description:
    A KBase module: COBRApyBasedFBA
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "git@github.com:kbaseapps/COBRApyBasedFBA.git"
    GIT_COMMIT_HASH = "30c47badca8a2a8808f0fd4e2373fc807a213565"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        #END_CONSTRUCTOR
        pass


    def run_fba_pipeline(self, ctx, params):
        """
        Run flux balance analysis and return ID of FBA object with results
        :param params: instance of type "RunFBAPipelineParams" -> structure:
           parameter "fbamodel_id" of type "fbamodel_id" (A string
           representing a FBAModel id.), parameter "fbamodel_workspace" of
           type "workspace_name" (A string representing a workspace name.),
           parameter "media_id" of type "media_id" (A string representing a
           Media id.), parameter "media_workspace" of type "workspace_name"
           (A string representing a workspace name.), parameter
           "target_reaction" of type "reaction_id" (A string representing a
           reaction id.), parameter "fba_output_id" of type "fba_id" (A
           string representing a FBA id.), parameter "workspace" of type
           "workspace_name" (A string representing a workspace name.),
           parameter "fva" of type "bool" (A binary boolean), parameter
           "minimize_flux" of type "bool" (A binary boolean), parameter
           "simulate_ko" of type "bool" (A binary boolean), parameter
           "all_reversible" of type "bool" (A binary boolean), parameter
           "feature_ko_list" of list of type "feature_id" (A string
           representing a feature id.), parameter "reaction_ko_list" of list
           of type "reaction_id" (A string representing a reaction id.),
           parameter "media_supplement_list" of list of type "compound_id" (A
           string representing a compound id.), parameter
           "objective_fraction" of Double, parameter "max_c_uptake" of
           Double, parameter "max_n_uptake" of Double, parameter
           "max_p_uptake" of Double, parameter "max_s_uptake" of Double,
           parameter "max_o_uptake" of Double, parameter "default_max_uptake"
           of Double
        :returns: instance of type "RunFBAPipelineResults" -> structure:
           parameter "new_fba_ref" of type "ws_fba_id" (The workspace ID for
           a FBA data object. @id ws KBaseFBA.FBA), parameter "objective" of
           Long, parameter "report_name" of String, parameter "report_ref" of
           type "ws_report_id" (The workspace ID for a Report object @id ws
           KBaseReport.Report)
        """
        # ctx is the context object
        # return variables are: results
        #BEGIN run_fba_pipeline
        
        kbase = cobrakbase.KBaseAPI(ctx['token'])
        fbamodel_json = kbase.get_object(params['fbamodel_id'], params['fbamodel_workspace'])
        fbamodel = cobrakbase.core.model.KBaseFBAModel(fbamodel_json)
        
        # Retrieve media
        media_json = kbase.get_object(params['media_id'], params['media_workspace'])
        media = cobrakbase.core.KBaseBiochemMedia(media_json)

        # TODO: add extra compounds to media with params['media_supplement_list']
        #       see what these params look like

        builder = KBaseFBAModelToCobraBuilder(fbamodel)
        model = builder.with_media(media).build()

        pipeline = FBAPipeline.fromKBaseParams(mock_params)
        result = pipeline.run(model, media)

        results = {'result': result}

        # TODO: save objects to ws
        #       Save FBA solution to KBase
        #       fba_output_id (id inside output result)
        #       make html report with report client 'sdk install <report client name>'
        
        #END run_fba_pipeline

        # At some point might do deeper type checking...
        if not isinstance(results, dict):
            raise ValueError('Method run_fba_pipeline return value ' +
                             'results is not type dict as required.')
        # return the results
        return [results]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
