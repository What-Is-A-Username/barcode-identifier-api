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

import os
import sys
import time
from typing import List
import requests
import platform
from xmltramp2 import xmltramp
from optparse import OptionParser
from Bio.SeqIO import SeqRecord

try:
    from urllib.parse import urlparse, urlencode
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError
    from urllib.request import __version__ as urllib_version
except ImportError:
    from urlparse import urlparse
    from urllib import urlencode
    from urllib2 import urlopen, Request, HTTPError
    from urllib2 import __version__ as urllib_version

try:
    unicode('')
except NameError:
    unicode = str

# Set interval for checking status
pollFreq = 3
# Output level
outputLevel = 1
# Debug level
debugLevel = 0
# Number of option arguments.
numOpts = len(sys.argv)

debugLevel = 0

base_url = 'https://www.ebi.ac.uk/Tools/services/rest/clustalo'

def printClustalODebugMessage(functionName, message, level):
    if (level <= debugLevel):
        print(u'[' + functionName + u'] ' + message, file=sys.stderr)

# User-agent for request (see RFC2616).
def getUserAgent():
    printClustalODebugMessage(u'getUserAgent', u'Begin', 11)
    # Agent string for urllib2 library.
    urllib_agent = u'Python-urllib/%s' % urllib_version
    clientRevision = '2023-02-01'
    # Prepend client specific agent string.
    try:
        pythonversion = platform.python_version()
        pythonsys = platform.system()
    except ValueError:
        pythonversion, pythonsys = "Unknown", "Unknown"

    # generate user_agent
    # example: 'EBI-Sample-Client/2022-09-13 12:15 (clustalo.py; Python 3.8.10; Linux) Python-urllib/3.8'
    user_agent = u'Barcode-Identifier-API-Server/%s (%s; Python %s; %s) %s' % (
        clientRevision, os.path.basename(__file__),
        pythonversion, pythonsys, urllib_agent)

    printClustalODebugMessage(u'getUserAgent', u'user_agent: ' + user_agent, 12)
    printClustalODebugMessage(u'getUserAgent', u'End', 11)
    return user_agent

def submitMultipleAlignmentAsync(sequence: str, run_id: str) -> str:    
    '''
    Submit a multiple alignment job
    '''
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

    print(params)

    requestUrl = f'{base_url}/run/'
    requestData = urlencode(params)

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
    printClustalODebugMessage(u'serviceRun', u'jobId: ' + jobId, 2)
    printClustalODebugMessage(u'serviceRun', u'End', 1)
    return jobId

def restRequest(url):
    printClustalODebugMessage(u'restRequest', u'Begin', 11)
    printClustalODebugMessage(u'restRequest', u'url: ' + url, 11)
    try:
        # Set the User-agent.
        user_agent = getUserAgent()
        http_headers = {u'User-Agent': user_agent}
        req = Request(url, None, http_headers)
        # Make the request (HTTP GET).
        reqH = urlopen(req)
        resp = reqH.read()
        contenttype = reqH.info()

        if (len(resp) > 0 and contenttype != u"image/png;charset=UTF-8"
                and contenttype != u"image/jpeg;charset=UTF-8"
                and contenttype != u"application/gzip;charset=UTF-8"):
            try:
                result = unicode(resp, u'utf-8')
            except UnicodeDecodeError:
                result = resp
        else:
            result = resp
        reqH.close()
    # Errors are indicated by HTTP status codes.
    except HTTPError as ex:
        result = requests.get(url).content
    printClustalODebugMessage(u'restRequest', u'End', 11)
    return result

def serviceGetStatus(jobId):
    '''
        Determine whether the job is completed by sending a REST request to ClustalO
    '''
    printClustalODebugMessage(u'serviceGetStatus', u'Begin', 1)
    requestUrl = base_url + u'/status/' + jobId
    result = restRequest(requestUrl)
    printClustalODebugMessage(u'serviceGetStatus', u'status: ' + result, 2)
    printClustalODebugMessage(u'serviceGetStatus', u'End', 1)
    print(result)
    return result == "FINISHED"

def serviceGetResult(jobId, type_):
    '''
        Get an individual result type
    '''
    printClustalODebugMessage(u'serviceGetResult', u'Begin', 1)
    printClustalODebugMessage(u'serviceGetResult', u'jobId: ' + jobId, 2)
    printClustalODebugMessage(u'serviceGetResult', u'type_: ' + type_, 2)
    requestUrl = base_url + u'/result/' + jobId + u'/' + type_
    result = restRequest(requestUrl)
    printClustalODebugMessage(u'serviceGetResult', u'End', 1)
    return result

# Get available result types for job
def serviceGetResultTypes(jobId):
    printClustalODebugMessage(u'serviceGetResultTypes', u'Begin', 1)
    printClustalODebugMessage(u'serviceGetResultTypes', u'jobId: ' + jobId, 2)
    requestUrl = base_url + u'/resulttypes/' + jobId
    printClustalODebugMessage(u'serviceGetResultTypes', u'requestUrl: ' + requestUrl, 2)
    xmlDoc = restRequest(requestUrl)
    doc = xmltramp.parse(xmlDoc)
    printClustalODebugMessage(u'serviceGetResultTypes', u'End', 1)
    return doc[u'type':]

def getMultipleAlignmentResult(job_id: str, run_id: str):
    '''
        Save results to local file system
    '''

    printClustalODebugMessage(u'getResult', u'Begin', 1)
    printClustalODebugMessage(u'getResult', u'jobId: ' + job_id, 1)
    if outputLevel > 1:
        print("Getting results for job %s" % job_id)

    # Get available result types
    resultTypes = serviceGetResultTypes(job_id)

    outfile = f'./runs/{run_id}/{job_id}'

    for resultType in resultTypes:
        # Derive the filename for the result
        filename = (outfile + u'.' + unicode(resultType[u'identifier']) +
                    u'.' + unicode(resultType[u'fileSuffix']))

        outformat_parm = str(None).split(',')
        for outformat_type in outformat_parm:
            outformat_type = outformat_type.replace(' ', '')

            if outformat_type == 'None':
                outformat_type = None

            if not outformat_type or outformat_type == unicode(resultType[u'identifier']):
                if outputLevel > 1:
                    print("Getting %s" % unicode(resultType[u'identifier']))
                # Get the result
                result = serviceGetResult(job_id, unicode(resultType[u'identifier']))
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
                if outputLevel > 0:
                    print("Creating result file: " + filename)
    printClustalODebugMessage(u'getResult', u'End', 1)

