# -*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser

from COBRApyBasedFBA.COBRApyBasedFBAImpl import COBRApyBasedFBA
from COBRApyBasedFBA.COBRApyBasedFBAServer import MethodContext
from COBRApyBasedFBA.authclient import KBaseAuth as _KBaseAuth

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.WorkspaceClient import Workspace


class COBRApyBasedFBATest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('COBRApyBasedFBA'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'COBRApyBasedFBA',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = COBRApyBasedFBA(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ContigFilter_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa
        cls.prepareTestData()

    @classmethod
    def prepareTestData(cls):
        """This function creates an assembly object for testing"""
        pass

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def test_your_method(self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        # 
        # params = {
        #     'fbamodel_id':           'Lactococcus_lactis_model',
        #     'fbamodel_workspace':    'chenry:narrative_1574200759205',
        #     'media_id':              'RichDefinedAnaerobic.media2',
        #     'media_workspace':       'chenry:narrative_1574200759205',
        #     'target_reaction':       '',  # TODO: What is this? of type "reaction_id" (A string representing a reaction id.)
        #                                     #       objective function to optimize
        #     'fba_output_id':         'result_fba_sol',  # TODO: What is this? of type "fba_id" (A string representing a FBA id.)
        #                                     # Name of object we are going to save (to put back in workspace)
        #     'workspace':             self.wsName,
        #     'fva':                   False,
        #     'minimize_flux':         False, # pfba
        #     'simulate_ko':           False,
        #     'all_reversible':        False,
        #     'feature_ko_list':       [], #['L192589', 'L182555'],     # TODO: what is a feature?
        #     'reaction_ko_list':      [], #['rxn05625_c0', 'rxn00545_c0'],
        #     'media_supplement_list': [],     # list of compound id's
        #     'objective_fraction':    1.,
        #     'max_c_uptake':          0.,
        #     'max_n_uptake':          0.,
        #     'max_p_uptake':          0.,
        #     'max_s_uptake':          0.,
        #     'max_o_uptake':          0.,
        #     'default_max_uptake':    0.
        # }

        params = {
            'fbamodel_id':           'Ec_core_flux1',
            'fbamodel_workspace':    'abrace05:narrative_1594056508275',
            'media_id':              'Carbon-D-Glucose',
            'media_workspace':       'abrace05:narrative_1594056508275',
            'target_reaction':       '',  # TODO: What is this? of type "reaction_id" (A string representing a reaction id.)
                                            #       objective function to optimize
            'fba_output_id':         'result_fba_sol',  # TODO: What is this? of type "fba_id" (A string representing a FBA id.)
                                            # Name of object we are going to save (to put back in workspace)
            'workspace':             self.wsName,
            'fva':                   False,
            'minimize_flux':         False, # pfba # TODO: these are ints 0, 1
            'simulate_ko':           False,
            'all_reversible':        False,
            'feature_ko_list':       [], #['L192589', 'L182555'],     # TODO: what is a feature?
            'reaction_ko_list':      [], #['rxn05625_c0', 'rxn00545_c0'],
            'media_supplement_list': [],     # list of compound id's
            'objective_fraction':    1.,
            'max_c_uptake':          0.,
            'max_n_uptake':          0.,
            'max_p_uptake':          0.,
            'max_s_uptake':          0.,
            'max_o_uptake':          0.,
            'default_max_uptake':    0.
        }

        # From kbase
        params = {'workspace': 'abrace05:narrative_1594056508275', 'fbamodel_id': '44773/6/1', 'media_id': '44773/2/1', 'target_reaction': 'bio1', 'fba_output_id': 'FBA_test_result', 'fva': 0, 'minimize_flux': 0, 'simulate_ko': 0, 'feature_ko_list': 'b0001', 'reaction_ko_list': '', 'media_supplement_list': '', 'max_c_uptake': None, 'max_n_uptake': None, 'max_p_uptake': None, 'max_s_uptake': None, 'max_o_uptake': None}
        self.serviceImpl.run_fba_pipeline(self.ctx, params)
