#!/usr/bin/python2.7

# ----------------------------------------------------------------------#
# Copyright (c) 2013, Kaushalya Amarasinghe. 2022, Stephen G. Gaffney
#
# > Source License <
# This file is part of ADTEx.
#
#    ADTEx is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ADTEx is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with ADTEx.  If not, see <http://www.gnu.org/licenses/>.
#
#
# -----------------------------------------------------------------------#

import os
import sys
import shlex
import shutil
import argparse
import datetime
import subprocess
from multiprocessing import Process

from getMeanCoverage import getMeanCoverage

# absolute script path
scriptPath = os.path.realpath(os.path.dirname(sys.argv[0]))


class Options:
    def __init__(self):
        self.parser = argparse.ArgumentParser("Aberration Detection in Tumour EXome")
        self.parser.add_argument(
            "-n",
            "--normal",
            help="Matched normal sample in BAM format or bed formatted coverage [REQUIRED], "
            "to generate bed formatted coverage please see documentation",
            dest="control",
        )
        self.parser.add_argument(
            "-t",
            "--tumor",
            help="Tumor sample in BAM format or bed format DOC coverage [REQUIRED], "
            "to generate bed formatted coverage please see documentation",
            dest="tumor",
        )
        self.parser.add_argument(
            "-b",
            "--bed",
            help="BED format of the targeted regions [REQUIRED]",
            dest="bed",
        )
        self.parser.add_argument(
            "-o",
            "--out",
            help="Output folder path name to store the output of analysis [REQUIRED]",
            action="store",
            dest="outFolder",
        )
        self.parser.add_argument(
            "--DOC",
            help="If specified, matched normal and tumor inputs  will be in BED "
            "formatted coverage [False]",
            action="store_true",
            dest="doc",
            default="False",
        )
        self.parser.add_argument(
            "--ploidy",
            help="Most common ploidy in the tumour sample [2]",
            dest="ploidy",
            action="store",
        )
        self.parser.add_argument(
            "--estimatePloidy",
            help="If provided, --baf must be specified to estimate base ploidy [FALSE]",
            dest="p_est",
            action="store_true",
            default="False",
        )
        self.parser.add_argument(
            "--minReadDepth",
            help="The threshold for minimum read depth for each exon [10]",
            action="store",
            dest="minReadDepth",
            default=10,
        )
        self.parser.add_argument(
            "-p",
            "--plot",
            help="Plots each chromosome with CNV estimates [False]",
            action="store_true",
            dest="plot",
            default="False",
        )
        self.parser.add_argument(
            "--baf",
            help="File containing B allele frequencies at heterozygous loci of "
            "the normal [optional]",
            dest="baf",
            action="store",
        )

        args = self.parser.parse_args()
        if args.control:
            self.control = args.control
        else:
            self.parser.error("Matched normal sample not supplied")
        if args.tumor:
            self.tumor = args.tumor
        else:
            self.parser.error("Tumor sample not supplied")
        if args.bed:
            self.bed = args.bed
        else:
            self.parser.error("Targeted regions file not supplied")
        if args.baf:
            self.baf = args.baf
            self.bafin = "True"
        else:
            self.bafin = "False"
            self.baf = "False"
        if args.outFolder:
            self.outFolder = str(args.outFolder)
        else:
            self.parser.error("Output folder path not supplied")
        if str(args.doc) == "True":
            # self.totReads = args.doc
            self.docInput = "True"
        else:
            self.docInput = "False"
        if args.ploidy:
            self.ploidy = args.ploidy
            self.p_est = "False"
            if str(args.p_est) == "True":
                print("Ploidy provided. Estimation step won't be extecuted")
            self.ploidyIn = "True"
        else:
            self.ploidyIn = "False"
        if str(args.p_est) == "True" and str(self.ploidyIn) == "False":
            if args.baf:
                self.p_est = "True"
            else:
                self.parser.error(
                    "--baf must be provided if base ploidy estimation is used"
                )
        elif str(self.ploidyIn) == "False":
            self.ploidy = 2
            self.p_est = "False"
        if args.minReadDepth:
            self.minRead = args.minReadDepth
        if args.plot:
            self.plot = args.plot


def splitBam(inF, outFolder, chroms):
    try:
        os.mkdir(outFolder + "/chr/")
    except:
        print("Folder exists")

    for c in chroms:
        outputfile = outFolder + "/chr/" + c + ".bam"
        subprocess.call(
            "samtools view -bh %s %s > %s" % (inF, c, outputfile), shell=True
        )


def splitBed(inF, outFolder):
    try:
        os.mkdir(outFolder + "/chr/")
    except:
        print("Folder exists")
    inF = open(inF, "r")
    check = "1"
    outputfile = outFolder + "/chr/" + str(check) + ".bed"
    outfile = open(outputfile, "w")

    for row in inF:
        cols = row.split()
        chr = cols[0]
        if chr != check:
            outfile.close()
            check = chr
            outputfile = outFolder + "/chr/" + str(check) + ".bed"
            outfile = open(outputfile, "w")
        outfile.write(row.rstrip() + "\n")

    outfile.close()


def getCoverage(outF, bedF, chroms):
    outFile = outF + "/coverage.txt"
    iOutFile = open(outFile, "wb")
    for c in chroms:
        inFile = outF + "/chr/" + c + ".bam"
        targets = bedF + "/chr/" + c + ".bed"
        args = shlex.split("coverageBed -b %s -d -a %s" % (inFile, targets))
        output = subprocess.Popen(args, stdout=subprocess.PIPE).communicate()[0]
        iOutFile.write(output)
    iOutFile.close()
    subprocess.call(
        "sort -V -k1,1 -k2,2n -k3,3n -k4,4n %s > %s" % (outFile, outFile + ".sorted"),
        shell=True,
    )


def getCoveragefromBAM(outF, bedF, inF, genome_path):
    outFile = outF + "/coverage.txt"
    targets = bedF
    inFile = inF

    # to remove duplicates
    # outbam = outF+"/noduplicates.bam"
    # subprocess.call("samtools view -F 0x400 %s > %s" %(inFile,outbam),shell=True)
    # subprocess.call("coverageBed -abam %s -d -b %s > %s" %(outbam,targets,outFile),shell=True)

    subprocess.call(
        "coverageBed -b %s -d -a %s -g %s -sorted > %s"
        % (inFile, targets, genome_path, outFile),
        shell=True,
    )  # SGG: add genome path
    # subprocess.call("sort -V -k1,1 -k2,2n -k3,3n -k4,4n %s > %s" %(outFile,outFile+".sorted"),shell=True)
    shutil.copyfile(
        outFile, outFile + ".sorted"
    )  # SGG: copy file instead of sorting. assume already sorted.


def sortFile(inF, fileN):
    inFile = inF + fileN
    subprocess.call(
        "sort -V -k1,1 -k2,2n -k3,3n -k4,4n %s > %s" % (inFile, inFile + ".sorted"),
        shell=True,
    )


def getTotReads(inF, folder):
    subprocess.call(
        "samtools view %s | wc -l > %s/tot_reads.txt" % (inF, folder), shell=True
    )


def analyseCNV(params, ratio_data, outF, chroms):
    print("Analysing CNV...")
    rScriptName = os.path.join(scriptPath, "cnv_analyse.R")
    rFunctionsPath = os.path.join(scriptPath, "RFunction.R")

    def runCNV(par, ratios, outLoc, chrom, rScr, rFunctions, p):
        args = shlex.split(
            "Rscript %s %s %s %s %s %s %s %s %s %s"
            % (
                rScr,
                rFunctions,
                ratios,
                p,
                par.minRead,
                outLoc,
                par.bafin,
                par.baf,
                par.plot,
                chrom,
            )
        )
        _ = subprocess.call(args)

    if str(params.p_est) == "False":
        args = shlex.split(
            "Rscript %s %s %s %s %s %s %s %s %s %s"
            % (
                rScriptName,
                rFunctionsPath,
                ratio_data,
                params.ploidy,
                params.minRead,
                outF,
                params.bafin,
                params.baf,
                params.plot,
                chroms,
            )
        )
        _ = subprocess.call(args)
        args = shlex.split(
            "mv %s %s"
            % (outF + "/temp/cnv.result" + str(params.ploidy), outF + "/cnv.result")
        )
        _ = subprocess.call(args)
    else:
        cnv2 = Process(
            target=runCNV,
            args=(params, ratio_data, outF, chroms, rScriptName, rFunctionsPath, 2),
        )
        cnv3 = Process(
            target=runCNV,
            args=(params, ratio_data, outF, chroms, rScriptName, rFunctionsPath, 3),
        )
        cnv4 = Process(
            target=runCNV,
            args=(params, ratio_data, outF, chroms, rScriptName, rFunctionsPath, 4),
        )
        cnv2.start()
        cnv3.start()
        cnv4.start()
        cnv2.join()
        cnv3.join()
        cnv4.join()


def segmentRatio(params, cCoverage, tCoverage, outF, chroms) -> None:
    rScriptName = os.path.join(scriptPath, "segment_ratio.R")
    rFunctionsPath = os.path.join(scriptPath, "RFunction.R")
    args = shlex.split(
        "Rscript %s %s %s %s %s %s %s %s %s %s"
        % (
            rScriptName,
            rFunctionsPath,
            cCoverage,
            tCoverage,
            params.minRead,
            outF,
            params.bafin,
            params.baf,
            params.ploidyIn,
            chroms,
        )
    )
    print(f"R subprocess:\n {' '.join(args)}")
    subprocess.call(args)


def zygosity(params, outF, chroms) -> None:
    print("Predicting Zygosity states")
    rScriptName = os.path.join(scriptPath, "zygosity.R")
    args = shlex.split(
        "Rscript %s %s %s %s %s %s"
        % (rScriptName, params.baf, outF, params.minRead, params.plot, chroms)
    )
    print(f"R subprocess:\n {' '.join(args)}")
    subprocess.call(args)


def getChroms(inF, outF) -> str:
    outFile = outF + "/targets.sorted"
    # subprocess.call("sort -V -k1,1 -k2,2n -k3,3n -k4,4n %s > %s" %(inF,outFile),shell=True)
    # shutil.copyfile(inF, outFile)  # SGG: copy file instead of sorting. assume already sorted.
    os.symlink(
        get_full_path_sgg(inF), outFile
    )  # SGG: symlink file instead of sorting or copying.
    infile = open(outFile)
    chr = []
    chr_prev = 0
    for line in infile:
        chr_current = line.split("\t")[0]  # SGG: don't need rstrip
        if chr_current != chr_prev and chr_current[0] != "G":
            chr_prev = chr_current
            chr.append(chr_current)
    chr = ",".join(chr)
    return chr


def get_full_path_sgg(file_path) -> str:
    if file_path.startswith("/"):
        return file_path
    wd = os.getcwd()
    return os.path.join(wd, file_path)


def main():
    subprocess.call("date", shell=True)
    options = Options()
    control = options.control
    tumor = options.tumor
    targets = options.bed
    outF = options.outFolder
    docInput = options.docInput
    bafIn = options.bafin

    print("Creating output folder")

    if outF[len(outF) - 1] == "/":
        outF = outF[: len(outF) - 1]
    if os.path.exists(outF):
        now_str = datetime.datetime.utcnow().strftime('%Y-%m-%d_%H%M%S')
        shutil.move(outF, f"{outF}_{now_str}")
    os.mkdir(outF)


    os.mkdir(outF + "/temp")
    os.mkdir(outF + "/temp/control")
    os.mkdir(outF + "/temp/tumor")

    chroms = getChroms(targets, outF)
    targets = outF + "/targets.sorted"

    if str(docInput) == "True":
        print("Generating mean coverage files...")
        # SGG: symlink instead of copy. skip sorting step
        # subprocess.call("cp %s %s" %(control,outF+"/temp/control/coverage.txt"),shell=True)
        # subprocess.call("cp %s %s" %(tumor,outF+"/temp/tumor/coverage.txt"),shell=True)
        os.symlink(
            get_full_path_sgg(control), outF + "/temp/control/coverage.txt.sorted"
        )  # SGG
        os.symlink(
            get_full_path_sgg(tumor), outF + "/temp/tumor/coverage.txt.sorted"
        )  # SGG
        # subprocess.call("cp %s %s" %(options.totReads,outF+"/totReads.txt"),shell=True)
        # ctrSort=Process(target= sortFile, args=(outF+"/temp/control","/coverage.txt"))
        # tmrSort=Process(target= sortFile, args=(outF+"/temp/tumor","/coverage.txt"))
        # ctrSort.start()
        # tmrSort.start()
        # ctrSort.join()
        # tmrSort.join()

    else:
        print("Creating coverage files")
        # ctrSplit = Process(target = splitBam, args=(control,outF+"/temp/control",chroms))
        # tmrSplit = Process(target = splitBam, args=(tumor,outF+"/temp/tumor",chroms))
        # tSplit = Process(target = splitBed, args=(targets,outF+"/temp",chroms))
        # #ctrTotal = Process(target = getTotReads, args=(control,outF+"/temp/control"))
        # #tTotal = Process(target = getTotReads, args=(tumor,outF+"/temp/tumor"))

        # ctrSplit.start()
        # tmrSplit.start()
        # tSplit.start()
        # #ctrTotal.start()
        # #tTotal.start()

        # ctrSplit.join()
        # tmrSplit.join()
        # tSplit.join()
        # #ctrTotal.join()
        # #tTotal.join()

        # ctrDOC = Process(target= getCoverage, args=(outF+"/temp/control",outF+"/temp",chroms))
        # tmrDOC = Process(target= getCoverage, args=(outF+"/temp/tumor",outF+"/temp",chroms))
        # ctrDOC.start()
        # tmrDOC.start()
        # ctrDOC.join()
        # tmrDOC.join()

        # Following codes are for generating coverage for the whole bam at once
        genome_path = os.path.join(outF, os.path.basename(control) + "_genome.txt")
        # SGG: create genome.txt file
        """samtools view -H {bam} | grep -P "@SQ\tSN:" | sed 's/@SQ\tSN://' | sed 's/\tLN:/\t/' > {genome_path}"""
        cmd1 = "samtools view -H {bam}".format(bam=control)
        cmd2 = 'grep -P "@SQ\tSN:"'
        cmd3 = "sed 's/@SQ\tSN://'"
        cmd4 = "sed 's/\tLN:/\t/'"

        with open(genome_path, "w") as genome_file:
            p1 = subprocess.Popen(shlex.split(cmd1), stdout=subprocess.PIPE)
            p2 = subprocess.Popen(
                shlex.split(cmd2), stdout=subprocess.PIPE, stdin=p1.stdout
            )
            p3 = subprocess.Popen(
                shlex.split(cmd3), stdout=subprocess.PIPE, stdin=p2.stdout
            )
            p4 = subprocess.Popen(
                shlex.split(cmd4), stdout=genome_file, stdin=p3.stdout
            )
        p4.communicate()

        ctrDOC = Process(
            target=getCoveragefromBAM,
            args=(outF + "/temp/control", targets, control, genome_path),
        )
        tmrDOC = Process(
            target=getCoveragefromBAM,
            args=(outF + "/temp/tumor", targets, tumor, genome_path),
        )
        ctrDOC.start()
        tmrDOC.start()
        ctrDOC.join()
        tmrDOC.join()

    ctrDOC = Process(
        target=getMeanCoverage,
        args=(outF + "/temp/control/coverage.txt.sorted", outF + "/control.coverage"),
    )
    tmrDOC = Process(
        target=getMeanCoverage,
        args=(outF + "/temp/tumor/coverage.txt.sorted", outF + "/tumor.coverage"),
    )
    ctrDOC.start()
    tmrDOC.start()
    ctrDOC.join()
    tmrDOC.join()

    subprocess.call("rm -rf %s" % (outF + "/temp"), shell=True)

    os.mkdir(outF + "/temp")

    segmentRatio(
        options, outF + "/control.coverage", outF + "/tumor.coverage", outF, chroms
    )

    chroms = open(outF + "/chrom").readline()

    print("running analyse CNV")  # SGG print
    analyseCNV(options, outF + "/ratio.data", outF, chroms)

    if str(options.p_est) == "True":
        print("Estimating base ploidy...")
        rScriptName = os.path.join(scriptPath, "extract_cnv.R")
        args = shlex.split("Rscript %s %s" % (rScriptName, outF + "/temp"))
        print(f"R subprocess:\n {' '.join(args)}")
        _ = subprocess.call(args)

        print("intersecting snp segments with cnv tables")  # SGG print
        cmd = "intersectBed -a %s -b %s -wb > %s" % (
                outF + "/temp/snp_segments",
                outF + "/temp/cnv2",
                outF + "/temp/cnv2_baf.txt",
            )
        print(f"R subprocess:\n {cmd}")
        subprocess.call(cmd, shell=True)

        cmd = "intersectBed -a %s -b %s -wb > %s" % (
            outF + "/temp/snp_segments",
            outF + "/temp/cnv3",
            outF + "/temp/cnv3_baf.txt",
        )
        print(f"R subprocess:\n {cmd}")
        subprocess.call(cmd, shell=True)

        cmd = "intersectBed -a %s -b %s -wb > %s" % (
                outF + "/temp/snp_segments",
                outF + "/temp/cnv4",
                outF + "/temp/cnv4_baf.txt",
            )
        print(f"R subprocess:\n {cmd}")
        subprocess.call(cmd, shell=True)

        rScriptName = os.path.join(scriptPath, "base_cnv.R")
        args = shlex.split("Rscript %s %s" % (rScriptName, outF + "/temp"))
        print(f"R subprocess:\n {' '.join(args)}")
        _ = subprocess.call(args)

        ploidy = open(outF + "/temp/ploidy").readline()

        print("working on zygosity")  # SGG print
        args = shlex.split(
            "mv %s %s" % (outF + "/temp/cnv.result" + str(ploidy), outF + "/cnv.result")
        )
        print(f"R subprocess:\n {' '.join(args)}")
        _ = subprocess.call(args)

    print("removing temp dir")  # SGG print
    subprocess.call("rm -rf %s" % (outF + "/temp"), shell=True)

    if str(options.plot) == "True":
        rScriptName = os.path.join(scriptPath, "plot_results.R")
        args = shlex.split("Rscript %s %s %s" % (rScriptName, outF, chroms))
        print(f"R subprocess:\n {' '.join(args)}")
        _ = subprocess.call(args)

    print("working on zygosity")  # SGG print
    if str(bafIn) == "True":
        subprocess.call("mkdir %s" % (outF + "/zygosity"), shell=True)
        zygosity(options, outF, chroms)
        if str(options.plot) == "True":
            subprocess.call(
                "ls -v %s > %s"
                % (outF + "/zygosity/*.png", outF + "/zygosity/filelist"),
                shell=True,
            )
            l = open(outF + "/zygosity/filelist").read().split("\n")
            l = l[0 : (len(l) - 1)]
            l.insert(0, "convert")
            l.append(outF + "/zygosity/zygosity_results.pdf")
            print(f"R subprocess:\n {' '.join(l)}")
            subprocess.call(l)
            subprocess.call("rm %s" % (outF + "/zygosity/filelist"), shell=True)

    # subprocess.call("rm %s %s %s" %(outF+"/chrom",outF+"/ratio.data",outF+"/targets.sorted"),shell=True)  # SGG hide

    subprocess.call("date", shell=True)


if __name__ == "__main__":
    main()
