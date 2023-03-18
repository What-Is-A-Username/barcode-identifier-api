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
# Clustal Omega (REST) web service Python client using xmltramp2.
#
# For further information see:
# https://www.ebi.ac.uk/Tools/webservices/
#
###############################################################################

from __future__ import print_function

import sys
import time
from xmltramp2 import xmltramp
from barcode_blastn.file_paths import get_data_run_path

from barcode_blastn.helper.embl_utils import printDebugMessage, getUserAgent, restRequest, serviceGetResult, serviceGetResultTypes

try:
    from urllib.parse import urlencode
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError
except ImportError:
    from urllib import urlencode
    from urllib2 import urlopen, Request, HTTPError

try:
    unicode('')
except NameError:
    unicode = str

# Set interval for checking status
pollFreq = 2
# Output level
outputLevel = 1
# Debug level
debugLevel = 0
# Number of option arguments.
numOpts = len(sys.argv)

debugLevel = 0

base_url = 'https://www.ebi.ac.uk/Tools/services/rest/clustalo'

# TODO: Make call sync instead of async
def submitMultipleAlignmentAsync(sequence: str, run_id: str) -> str:  
    '''
    Submit a multiple alignment job
    '''

    print(f'Starting multiple alignment job for run {run_id}')

    params = {
        'email': 'william.huang1212@gmail.com',
        'title': f'run-{run_id}',
        'sequence': sequence, 
        'guidetreeout': 'true', 
        'addformats': 'false', 
        'dismatout': 'false', 
        'dealign': 'false', 
        'mbed': 'true', 
        'mbediteration': 'true', 
        'iterations': 0, 
        'gtiterations': -1, 
        'hmmiterations': -1, 
        'outfmt': 'clustal_num',
        'order': 'aligned'
    }
    
    requestUrl = f'{base_url}/run/'
    requestData = urlencode(params)

    try:
        # Set the HTTP User-agent.
        user_agent = getUserAgent()
        http_headers = {u'User-Agent': user_agent}
        req = Request(requestUrl, None, http_headers)
        # Make the submission (HTTP POST).
        print("Sending job to server.")
        reqH = urlopen(req, requestData.encode(encoding=u'utf_8', errors=u'strict'))
        jobId = unicode(reqH.read(), u'utf-8')
        print(f"Received job id {jobId}")
        reqH.close()
    except HTTPError as ex:
        print(xmltramp.parse(unicode(ex.read(), u'utf-8'))[0][0])
        quit()
    printDebugMessage(u'serviceRun', u'jobId: ' + jobId, 2)
    printDebugMessage(u'serviceRun', u'End', 1)
    return jobId

def serviceGetStatus(jobId):
    '''
        Determine whether the job is completed by sending a REST request to ClustalO
    '''
    printDebugMessage(u'serviceGetStatus', u'Begin', 1)
    requestUrl = base_url + u'/status/' + jobId
    result = restRequest(requestUrl)
    printDebugMessage(u'serviceGetStatus', u'status: ' + result, 2)
    printDebugMessage(u'serviceGetStatus', u'End', 1)
    print(result)
    return result

# Client-side poll
def clientPoll(jobId):
    time.sleep(1)
    printDebugMessage(u'clientPoll', u'Begin', 1)
    result = u'PENDING'
    while result == u'RUNNING' or result == u'PENDING':
        result = serviceGetStatus(jobId)
        print("polling: ", result)
        if result == u'RUNNING' or result == u'PENDING':
            time.sleep(pollFreq)
    printDebugMessage(u'clientPoll', u'End', 1)

def getMultipleAlignmentResult(job_id: str, run_id: str):
    '''
        Save results to local file system
    '''

    printDebugMessage(u'getResult', u'Begin', 1)
    printDebugMessage(u'getResult', u'jobId: ' + job_id, 1)
    print("Getting ClustalO results for job %s" % job_id)

    # Poll
    clientPoll(job_id)

    print("Polling completed %s" % job_id)

    # Get available result types
    resultTypes = serviceGetResultTypes(base_url, job_id)

    outfile = get_data_run_path(run_id=run_id) + '/' + job_id

    for resultType in resultTypes:
        # Derive the filename for the result
        filename = (outfile + u'.' + unicode(resultType[u'identifier']) +
                    u'.' + unicode(resultType[u'fileSuffix']))

        # TODO: See if the number of output files transferred can be reduced
        outformat_parm = str(None).split(',')
        for outformat_type in outformat_parm:
            outformat_type = outformat_type.replace(' ', '')

            if outformat_type == 'None':
                outformat_type = None

            if not outformat_type or outformat_type == unicode(resultType[u'identifier']):
                if outputLevel > 1:
                    print("Getting %s" % unicode(resultType[u'identifier']))
                # Get the result
                result = serviceGetResult(job_id, base_url, unicode(resultType[u'identifier']))
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

