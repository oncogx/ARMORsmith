import xml.etree.ElementTree as ET
import csv, sys, getopt
import pandas as pd
import re
import time
import os

import datetime
import numpy as np
from pathlib import Path

#Program to parse Foundation Medicine XML formatted data into csv format for PREDICT Server
#inputdir : the directory where the xmls are located
#outputdir : where to output the tables
#pathtopredict : location of the PREDICT - MRN table
#Note that error logs are written to file in parser directory
inputdir = '/home/scp010/ARMOR/data/foundmed/'
outputdir = '/home/scp010/ARMOR/data/test_processed/'

#pathtohistory = './parsehistory.txt'
overwrite = 0
currtime = time.strftime("%d-%m-%Y_%H-%M-%S")
HISTORY = []

def main(argv):
    global pathtohistory
    global HISTORY
    global inputdir
    global outputdir
    global overwrite
    try:
        opts, args = getopt.getopt(argv, "i:o:Ch", ["inputdir=","outputdir=","Clear","help"])
    except getopt.GetoptError:
        print('xmlparser.py -C <create new table> -i <inputdir> -o <outputdir>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('xmlparser.py -C <clear tables and create new table> -i <inputdir> -o <outputdir>')
            sys.exit()
        elif opt in ("-i", "--inputdir"):
            inputdir = str(arg)
        elif opt in ("-o", "--outputdir"):
            outputdir = str(arg)
        elif opt in ("-C", "--Clear"):
            overwrite = 1
            overwritecdwr = 1

    with open(pathtohistory) as fil:
        for row in csv.reader(fil, delimiter = '\t'):
            HISTORY.append(row[0].strip())
    for file in os.listdir(inputdir):
        if file.endswith(".xml"):

            line = str(os.path.join(inputdir, file))
            try:
                parseXML(line)
            except Exception as ex:
                errfilename = outputdir+'Error-File'+currtime
                nonparsedfilename = 'Not-Parsed'+currtime
                nonparsedfile = open(str(nonparsedfilename), 'a+')
                errfile = open(str(errfilename), 'a+')
                nonparsedfile.write(line)
                nonparsedfile.close()
                errfile.write(line)
                errfile.write('\n')
                errfile.write(str(ex))
                errfile.write('\n')
                errfile.close()


def parseXML(file):
    global hasRan
    global currtime
    global outputdir
    global pathtohistory
    global HISTORY
    global overwrite
    meta = pd.DataFrame(columns = ('MRN', 'TRF', 'Short-Variants', 'Rearrangements', 'Copy-Number-Alterations',
                                   'Non-Human-Content', 'First-Name', 'Last-Name', 'Gender', 'DOB', 'Record-Date', 'Collection-Date', 'Disease', 'Disease-Ontology',
                                   'Purity-Assessment', 'Specimen-Site', 'Quality', 'Test-Type','Tumor-Mutation-Burden', 'Tumor-Mutation-Burden-Score','Microsatellite-Instability'))
    mutations = pd.DataFrame(columns = ('MRN', 'TRF','Category', 'Status', 'Gene','Protein-Effect', 'Functional-Effect', 'Allele-Fraction',
                                        'Equivocal', 'Subclonal', 'Percent-Reads', 'CDS-Effect', 'Strand','Position', 'Copy-Number', 'Number-of-Exons',
                                        'Ratio', 'Description', 'inFrame', 'Other-Gene', 'pos1', 'pos2','supporting-read-pairs', 'organism', 'reads-per-million'))
    xmlns = {"rr":"http://integration.foundationmedicine.com/reporting"}
    tree = ET.parse(file.strip())
    root = tree.getroot()
    variantreport = root.find('.//{http://foundationmedicine.com/compbio/variant-report-external}variant-report')
    firstname = root.findtext('.//FirstName')
    lastname = root.findtext('.//LastName')
    mrn2 = root.findtext('.//MRN')
    trf = root.findtext('.//ReportId')
    if trf in HISTORY:
        return
    recDate = root.findtext('.//ReceivedDate')
    collDate = root.findtext('.//CollDate')
    DOB = root.findtext('.//DOB')
    disease = variantreport.get('disease')
    diseaseontology = variantreport.get('disease-ontology')
    gender = variantreport.get('gender')
    purityassessment = variantreport.get('purity-assessment')
    testtype = variantreport.get('test-type')
    specsite = variantreport.get('tissue-of-origin')
    quality = variantreport.find(".//{http://foundationmedicine.com/compbio/variant-report-external}quality-control").get('status', 'NaN')
    if variantreport.find(".//{http://foundationmedicine.com/compbio/variant-report-external}tumor-mutation-burden") is not None:
        tmbscore = variantreport.find(".//{http://foundationmedicine.com/compbio/variant-report-external}tumor-mutation-burden").get('score', 'NaN')
        tmbstatus = variantreport.find(".//{http://foundationmedicine.com/compbio/variant-report-external}tumor-mutation-burden").get('status', 'NaN')
    else:
        tmbscore = "NaN"
        tmbstatus = "NaN"
    if variantreport.find(".//{http://foundationmedicine.com/compbio/variant-report-external}microsatellite-instability") is not None:
        mss = variantreport.find(".//{http://foundationmedicine.com/compbio/variant-report-external}microsatellite-instability").get('status', 'NaN')
    else:
        mss = 'NaN'
    numsvariants = len(variantreport.findall(".//{http://foundationmedicine.com/compbio/variant-report-external}short-variant"))
    for svariant in variantreport.findall(".//{http://foundationmedicine.com/compbio/variant-report-external}short-variant"):
        mutations = mutations.append({ 'MRN': mrn2, 'TRF': trf,'Category': 'shortvariant',
                           'Status': svariant.get('status', 'NaN'),
                           'Gene': svariant.get('gene','NaN'),
                           'Protein-Effect':svariant.get('protein-effect','NaN'),
                           'Functional-Effect':svariant.get('functional-effect','NaN'),
                           'Allele-Fraction':svariant.get('allele-fraction','NaN'),
                           'Equivocal':svariant.get('equivocal','NaN'),
                           'Subclonal':svariant.get('subclonal','NaN'),
                           'Percent-Reads':svariant.get('percent-reads','NaN'),
                           'CDS-Effect':svariant.get('cds-effect','NaN'),
                           'Strand':svariant.get('strand','NaN'),
                           'Position':svariant.get('position','NaN')}, ignore_index = True)
    numrearrangements = len(variantreport.findall(".//{http://foundationmedicine.com/compbio/variant-report-external}rearrangement"))
    for rearrangement in variantreport.findall(".//{http://foundationmedicine.com/compbio/variant-report-external}rearrangement"):
        mutations = mutations.append({ 'MRN': mrn2, 'TRF': trf,'Category' : 'rearrangement',
                                       'Status': rearrangement.get('status', 'NaN'),
                                       'Allele-Fraction' : rearrangement.get('allele-fraction', 'NaN'),
                                       'Description' : rearrangement.get('description', 'NaN'),
                                       'Equivocal' : rearrangement.get('equivocal', 'NaN'),
                                       'inFrame' : rearrangement.get('in-frame', 'NaN'),
                                       'Other-Gene' : rearrangement.get('other-gene', 'NaN'),
                                       'Percent-Reads' : rearrangement.get('percent-reads','NaN'),
                                       'pos1' : rearrangement.get('pos1','NaN'),
                                       'pos2':rearrangement.get('pos2','NaN'),
                                       'Gene' :rearrangement.get('targeted-gene','NaN'),
                                       'supporting-read-pairs' : rearrangement.get('supporting-read-pairs','NaN'),
                                       'Functional-Effect': rearrangement.get('type','NaN')}, ignore_index = True)
    numcopynum = len(variantreport.findall(".//{http://foundationmedicine.com/compbio/variant-report-external}copy-number-alteration"))

    for copy in variantreport.findall(".//{http://foundationmedicine.com/compbio/variant-report-external}copy-number-alteration"):
        mutations = mutations.append({ 'MRN' : mrn2, 'TRF': trf,'Category' : 'copynumber',
                                       'Status': copy.get('status', 'NaN'),
                                       'Gene' : copy.get('gene','NaN'),
                                       'Equivocal':copy.get('equivocal','NaN'),
                                       'Number-of-Exons' : copy.get('number-of-exons','NaN'),
                                       'Position':copy.get('position','NaN'),
                                       'Ratio':copy.get('ratio','NaN'),
                                       'Functional-Effect':copy.get('type','NaN'),
                                       'Copy-Number':copy.get('copy-number','NaN')}, ignore_index = True)

    numnonhuman = len(variantreport.findall(".//{http://foundationmedicine.com/compbio/variant-report-external}non-human"))
    for nonhuman in variantreport.findall(".//{http://foundationmedicine.com/compbio/variant-report-external}non-human"):
            mutations = mutations.append({ 'MRN' : mrn2, 'TRF': trf,'Category' : 'nonhuman',
                                           'Status' : nonhuman.get('status', 'NaN'),
                                           'organism' : nonhuman.get('organism', 'NaN'),
                                           'reads-per-million' : nonhuman.get('reads-per-million', 'NaN')}, ignore_index = True)
    meta = meta.append({'MRN' : mrn2, 'TRF' : trf, 'Short-Variants' : numsvariants, 'Rearrangements' :numrearrangements , 'Copy-Number-Alterations' :numcopynum,
                                   'Non-Human-Content' :numnonhuman, 'First-Name' :firstname, 'Last-Name' :lastname, 'Gender' :gender , 'DOB' : DOB,
                        'Record-Date' :recDate, 'Collection-Date' :collDate, 'Disease' :disease, 'Disease-Ontology' :diseaseontology,
                                   'Purity-Assessment' :purityassessment, 'Specimen-Site' :specsite, 'Quality' :quality , 'Test-Type' :testtype,'Tumor-Mutation-Burden': tmbstatus, 'Tumor-Mutation-Burden-Score' : tmbscore,'Microsatellite-Instability' : mss}, ignore_index = True)

    summaryname = '' + outputdir + 'Summary'
    summarynamedate = '' + summaryname + currtime
    variantname = '' + outputdir + 'Variants'
    variantnamedate = '' + variantname + currtime
    if overwrite == 1:
        meta.to_csv(str(summaryname), sep = '\t', index = False, quoting = csv.QUOTE_ALL)
        mutations.to_csv(str(variantname), sep = '\t', index = False, quoting = csv.QUOTE_ALL)
        meta.to_csv(str(summarynamedate), sep = '\t', index = False,quoting = csv.QUOTE_ALL)
        mutations.to_csv(str(variantnamedate), sep = '\t', index = False,quoting = csv.QUOTE_ALL)
        overwrite = 0
    else:
        meta.to_csv(str(summaryname), mode = 'a+', header = False, sep = '\t', index = False, quoting = csv.QUOTE_ALL)
        mutations.to_csv(str(variantname), mode = 'a+', header = False, sep = '\t', index = False,quoting = csv.QUOTE_ALL)
        meta.to_csv(str(summarynamedate), mode = 'a+', header = False, sep = '\t', index = False,quoting = csv.QUOTE_ALL)
        mutations.to_csv(str(variantnamedate), mode = 'a+', header = False, sep = '\t', index = False,quoting = csv.QUOTE_ALL)

    historyfile = open(pathtohistory, 'a+')
    HISTORY.append(trf)
    historyfile.write("%s\n" % trf)
    historyfile.close()

main(sys.argv[1:])
