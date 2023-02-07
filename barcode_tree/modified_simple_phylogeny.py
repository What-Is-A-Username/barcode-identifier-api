#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#
# Copyright 2012-2022 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Python Client Automatically generated with:
# https://github.com/ebi-wp/webservice-clients-generator
#
# Simple Phylogeny (REST) web service Python client using xmltramp2.
#
# For further information see:
# https://www.ebi.ac.uk/Tools/webservices/
#
###############################################################################

from __future__ import print_function

import os
# TODO: see if another xml library can be used to replace xmltramp
from xmltramp2 import xmltramp

from barcode_tree.embl_utils import getUserAgent, printDebugMessage, restRequest, serviceGetResult, serviceGetResultTypes

try:
    from urllib.parse import urlencode
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError
# TODO: fix import errors
except ImportError:
    from urllib import urlencode
    from urllib2 import urlopen, Request, HTTPError

# allow unicode(str) to be used in python 3
try:
    unicode('')
except NameError:
    unicode = str

# Base URL for service
baseUrl = u'https://www.ebi.ac.uk/Tools/services/rest/simple_phylogeny'
version = u'2022-09-13 12:15'

def submitSimplePhylogenyAsync(run_id: str):
    # python simple_phylogeny.py --tree=phylip --clustering=Neighbour-joining runs/0ebab7ef-f7a8-456a-923c-a57b6e3e47ba/clustalo-R20230207-024801-0870-6020700-p1m.aln-clustal_num.clustal_num --email william.huang1212@gmail.com --asyncjob
    run_folder = os.path.abspath(f'./runs/{run_id}')
    run_contents = os.listdir(run_folder)
    try:
        alignment_file = [r for r in run_contents if r.endswith('.clustal_num')][0]
    except IndexError:
        print(f"Critical error: Could not find alignment file with extension .clustal_num for run {run_id}")
        raise FileNotFoundError()

    with open(f'{run_folder}/{alignment_file}') as alignment_handler:
        alignment_sequence = alignment_handler.read()
    params = { 
        'email': 'william.huang1212@gmail.com',
        'title': f'tree-{run_id}',
        'sequence': alignment_sequence,
        'tree': 'phylip', 
        'kimura': 'false', 
        'tossgaps': 'false', 
        'clustering': 'Neighbour-joining',
        'pim': 'false'
    }

    requestUrl = baseUrl + u'/run/'
    printDebugMessage(u'serviceRun', u'requestUrl: ' + requestUrl, 2)

    # Get the data for the other options
    requestData = urlencode(params)

    printDebugMessage(u'serviceRun', u'requestData: ' + requestData, 2)
    # Errors are indicated by HTTP status codes.
    try:
        # Set the HTTP User-agent.
        user_agent = getUserAgent()
        http_headers = {u'User-Agent': user_agent}
        req = Request(requestUrl, None, http_headers)
        # Make the submission (HTTP POST).
        reqH = urlopen(req, requestData.encode(encoding=u'utf_8', errors=u'strict'))
        jobId = unicode(reqH.read(), u'utf-8')
        reqH.close()
    except HTTPError as ex:
        print(xmltramp.parse(unicode(ex.read(), u'utf-8'))[0][0])
        quit()
    printDebugMessage(u'serviceRun', u'jobId: ' + jobId, 2)
    printDebugMessage(u'serviceRun', u'End', 1)
    print("Tree submission successful. Job Id is ", jobId)
    return jobId

def checkSimplePhylogenyStatus(jobId):
    '''
    Check the status of phylogeny construction job. Return true if finished.
    '''
    printDebugMessage(u'checkSimplePhylogenyStatus', u'Begin', 1)
    printDebugMessage(u'checkSimplePhylogenyStatus', u'jobId: ' + jobId, 2)
    requestUrl = baseUrl + u'/status/' + jobId
    printDebugMessage(u'checkSimplePhylogenyStatus', u'requestUrl: ' + requestUrl, 2)
    status = restRequest(requestUrl)
    printDebugMessage(u'checkSimplePhylogenyStatus', u'status: ' + status, 2)
    printDebugMessage(u'checkSimplePhylogenyStatus', u'End', 1)
    return status == 'FINISHED'

def getSimplePhylogenyOutput(jobId, run_id):
    printDebugMessage(u'getResult', u'Begin', 1)
    printDebugMessage(u'getResult', u'jobId: ' + jobId, 1)
    print("Getting Simple Phylogeny results for job %s" % jobId)

    # Get available result types
    resultTypes = serviceGetResultTypes(base_url=baseUrl, jobId=jobId)
    outfile = f'./runs/{run_id}/{jobId}'

    for resultType in resultTypes:
        # Derive the filename for the result
        filename = (outfile + u'.' + unicode(resultType[u'identifier']) +
                    u'.' + unicode(resultType[u'fileSuffix']))
        # Write a result file

        # TODO: See if the number of output files transferred can be reduced
        outformat_parm = str(None).split(',')
        for outformat_type in outformat_parm:
            outformat_type = outformat_type.replace(' ', '')

            if outformat_type == 'None':
                outformat_type = None

            if not outformat_type or outformat_type == unicode(resultType[u'identifier']):
                print("Getting %s" % unicode(resultType[u'identifier']))
                # Get the result
                result = serviceGetResult(jobId, baseUrl, unicode(resultType[u'identifier']))
                if (unicode(resultType[u'mediaType']) == u"image/png"
                        or unicode(resultType[u'mediaType']) == u"image/jpeg"
                        or unicode(resultType[u'mediaType']) == u"application/gzip"):
                    fmode = 'wb'
                else:
                    fmode = 'w'

                try:
                    fh = open(filename, fmode)
                    fh.write(result)
                    fh.close()
                except TypeError:
                    fh.close()
                    fh = open(filename, "wb")
                    fh.write(result)
                    fh.close()
                print("Creating result file: " + filename)
    printDebugMessage(u'getResult', u'End', 1)

def readTreeFromFile(run_id: str) -> str:
    '''
        Return the tree as a string by reading it from file
    '''
    run_folder = os.path.abspath(f'./runs/{run_id}')
    run_files = os.listdir(run_folder)
    try:
        phylip_file = [r for r in run_files if r.endswith(".tree.ph")][0]
    except IndexError:
        print("Warning: Tried to find local copy of tree PHYLIP file, but no file was found.")
        return ''
    else:
        with open(f'{run_folder}/{phylip_file}') as tree_file:
            tree_text = tree_file.read()
        return tree_text
