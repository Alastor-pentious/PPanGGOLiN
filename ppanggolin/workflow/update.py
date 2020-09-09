#!/usr/bin/env python3
#coding:utf-8

#default libraries
import os
import time
import argparse
import logging
from copy import deepcopy
#local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import mkFilename
from ppanggolin.annotate import annotatePangenome, readAnnotations, getGeneSequencesFromFastas
from ppanggolin.cluster import clustering, readClustering
from ppanggolin.graph import computeNeighborsGraph
from ppanggolin.nem.rarefaction import makeRarefactionCurve
from ppanggolin.nem.partition import partition
from ppanggolin.formats import writePangenome, getOrganismsNames, ErasePangenome, readPangenome
from ppanggolin.figures import drawTilePlot, drawUCurve
from ppanggolin.info import printInfo
from ppanggolin.RGP.genomicIsland import predictRGP
from ppanggolin.RGP.spot import predictHotspots
# a workflow to add genomes to an existing pangenome, and recompute said pangenome

def checkInputNames(currentNames, input):
    """Checks that no given input genome name conflicts with current genomes"""
    f = open(input,"r")
    conflictingNames = set()
    for line in f:
        name = [el.strip() for el in line.split("\t")][0]
        if name in currentNames:
            conflictingNames.add(name)
    if len(conflictingNames) != 0:
        raise Exception(f"Some of the given genomes have identical names to those that are already in the pangenome : '{' '.join(conflictingNames)}'")

def launch(args):
    #Hello courageous developper.
    #As for the futur, should we be capable of update the clustering step (instead of redoing everything) things can be ponctually updated rather than entirely recomputed, which can make an update much quicker than what exists now.
    #If you are here, it might be because you have been assigned this task. Good luck in your quest.

    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    if pangenome.parameters["cluster"]["read_clustering_from_file"]:
        raise Exception("This does not work for pangenomes whose clustering was provided previously. However, it is possible and was not been implemented as we did not have the use for it. If you're interested you can make yourself known on the PPanGGOLiN github page.")
    #checking that there is no conflict between given name and former names
    checkInputNames(getOrganismsNames(pangenome), args.fasta if args.fasta is not None else args.anno)

    #check conflicts between the way genomes are passed now, and the way they were passed before, and warn the user

    #annotate the new genomes and add them to the pangenome

    start_anno = time.time()

    if args.fasta is not None:
        
        #to cope with old pangenomes that did not have this option
        contig_filter = pangenome.parameters["annotation"].get("contig_filter")
        if contig_filter is None:
            contig_filter = 0
        annotatePangenome(pangenome, args.fasta, args.tmpdir, args.cpu, translation_table=pangenome.parameters["annotation"]["translation_table"], kingdom=pangenome.parameters["annotation"]["kingdom"], norna= not pangenome.parameters["annotation"]["annotate_RNA"], overlap = pangenome.parameters["annotation"]["remove_Overlapping_CDS"], contig_filter = contig_filter)

    if args.anno is not None:
        read_pseudo = pangenome.parameters["annotation"].get("use_pseudo", False)
        readAnnotations(pangenome, args.anno, args.cpu, pseudo=read_pseudo)

    annotime = time.time() - start_anno

    formerParameters = deepcopy(pangenome.parameters)
    formerStatuses = deepcopy(pangenome.status)

    readPangenome(pangenome, annotation = True, geneSequences=True)#read all of the old annotations and gene sequences

    start_writing = time.time()
    writePangenome(pangenome, pangenome.file, force=True)
    writing_time = time.time() - start_writing

    #now we have all our genomes annotated, but some are loaded and some are not. 
    # atm we cannot cherry pick genomes so We'll unload what's been computed so far and reload everything with the formerly used parameters
    #Erasing formerly computed pangenome
    ErasePangenome(pangenome, geneFamilies = True)
    
    #rerunning each workflow that was formerly computed with the parameters that were used.
    force = True
    if formerStatuses["genesClustered"] != "No":
        start_clust = time.time()
        clustering(pangenome, args.tmpdir, cpu=args.cpu, defrag=formerParameters["cluster"]["defragmentation"], code = formerParameters["cluster"]["translation_table"], coverage=pangenome.parameters["cluster"]["coverage"], identity=pangenome.parameters["cluster"]["identity"], force = force)
        clust_time = time.time() - start_clust
        if formerStatuses["neighborsGraph"] != "No":
            start_graph = time.time()
            computeNeighborsGraph(pangenome, remove_copy_number = formerParameters["graph"]["removed_high_copy_number_families"])
            graph_time = time.time() - start_graph
            if formerStatuses["partitionned"] != "No":
                start_part = time.time()
                if not formerParameters["partition"]["computed_K"]:
                    K = formerParameters["partition"]["K"]
                else:
                    K = -1
                partition(pangenome, tmpdir = args.tmpdir, cpu = args.cpu, K=K, beta = formerParameters["partition"]["beta"], sm_degree=formerParameters["partition"]["max_node_degree_for_smoothing"])
                part_time = time.time() - start_part

        start_writing = time.time()
        writePangenome(pangenome, pangenome.file, args.force)
        writing_time = writing_time + time.time() - start_writing

        if formerStatuses["predictedRGP"] != "No":
            start_regions = time.time()
            predictRGP(pangenome, persistent_penalty=formerParameters["RGP"]["persistent_penalty"],variable_gain=formerParameters["RGP"]["variable_gain"],min_length=formerParameters["RGP"]["min_length"], min_score=formerParameters["RGP"]["min_score"], dup_margin=formerParameters["RGP"]["dup_margin"], cpu = args.cpu )
            regions_time = time.time() - start_regions
            if formerStatuses["spots"] != "No":
                start_spots = time.time()
                predictHotspots(pangenome, None, cpu = args.cpu, overlapping_match=formerParameters["spots"]["overlapping_match"], set_size=formerParameters["spots"]["set_size"], exact_match=formerParameters["spots"]["exact_match"])
                spot_time = time.time() - start_spots

            start_writing = time.time()
            writePangenome(pangenome, pangenome.file, False)
            writing_time = writing_time + time.time() - start_writing

    logging.getLogger().info(f"Annotation : {round(annotime,2)} seconds")
    if formerStatuses["genesClustered"] != "No":
        logging.getLogger().info(f"Clustering : {round(clust_time,2)} seconds")
    if formerStatuses["neighborsGraph"] != "No":
        logging.getLogger().info(f"Building the graph : {round(graph_time,2)} seconds")
    if formerStatuses["partitionned"] != "No":
        logging.getLogger().info(f"Partitionning the pangenome : {round(part_time,2)} seconds")
    if formerStatuses["predictedRGP"] != "No":
        logging.getLogger().info(f"Predicting RGP : {round(regions_time,2)} seconds")
    if formerStatuses["spots"] != "No":
        logging.getLogger().info(f"Gathering RGP into spots : {round(spot_time,2)} seconds")
    logging.getLogger().info(f"Writing the pangenome data in HDF5 : {round(writing_time,2)} seconds")

    printInfo(pangenome.file, content = True)


def updateSubparser(subparser):
    parser = subparser.add_parser("update", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group(title = "Required arguments", description = "One of the following arguments is required :")
    required.add_argument('--fasta',  required=False, type=str, help="A tab-separated file listing the organism names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). One line per organism. This option can be used alone with the pangenome file.")
    required.add_argument('--anno', required=False, type=str, help="A tab-separated file listing the organism names, and the gff filepath of its annotations (the gffs can be compressed). One line per organism. This option can be used alone with the pangenome file IF the fasta sequences are in the gff files, otherwise --fasta needs to be used.")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file. This MUST be provided.")
     #required.add_argument("--clusters",required=False, type=str, help = "a tab-separated file listing the cluster names, the gene IDs, and optionnally whether they are a fragment or not.")
    return parser
