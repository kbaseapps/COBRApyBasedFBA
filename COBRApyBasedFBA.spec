/*
A KBase module: COBRApyBasedFBA
*/

module COBRApyBasedFBA {
    /*
        A binary boolean
    */
    typedef int bool;
    /*
        A string representing a Media id.
    */
    typedef string media_id;
    /*
        A string representing a FBAModel id.
    */
    typedef string fbamodel_id;
    /*
        A string representing a FBA id.
    */
    typedef string fba_id;
    /*
        A string representing a reaction id.
    */
    typedef string reaction_id;
    /*
        A string representing a feature id.
    */
    typedef string feature_id;
    /*
        A string representing a compound id.
    */
    typedef string compound_id;
    /*
        A string representing a workspace name.
    */
    typedef string workspace_name;
    /* 
        The workspace ID for a FBAModel data object.
        @id ws KBaseFBA.FBAModel
    */
    typedef string ws_fbamodel_id;
    /* 
        The workspace ID for a FBA data object.
        @id ws KBaseFBA.FBA
    */
    typedef string ws_fba_id;
    /* 
        The workspace ID for a Report object
        @id ws KBaseReport.Report
    */
    typedef string ws_report_id;

    typedef structure {
        reaction_id custom_reaction_id;
        float custom_lb;
        float custom_ub;
    } CustomBounds;
    
    typedef structure {
        fbamodel_id fbamodel_id;
        workspace_name fbamodel_workspace;
        media_id media_id;
        workspace_name media_workspace;
        reaction_id target_reaction;
        fba_id fba_output_id;
        workspace_name workspace;
        string solver;
        
        bool minimize_objective;
        string fba_type;
        string fva_type;
        bool simulate_ko;
        bool all_reversible;

        float fraction_of_optimum_pfba;
        float fraction_of_optimum_fva;
        
        list<feature_id> feature_ko_list;
        list<reaction_id> reaction_ko_list;
        list<compound_id> media_supplement_list;
        list<CustomBounds> custom_bound_list;

        float max_c_uptake;
        float max_n_uptake;
        float max_p_uptake;
        float max_s_uptake;
        float max_o_uptake;
        float default_max_uptake;
    } RunFBAPipelineParams;
    
    typedef structure {
        ws_fba_id new_fba_ref;
        int objective;
        string report_name;
        ws_report_id report_ref;
    } RunFBAPipelineResults;
    /*
        Run flux balance analysis and return ID of FBA object with results 
    */
    funcdef run_fba_pipeline(RunFBAPipelineParams params) returns (RunFBAPipelineResults results) authentication required;

};
