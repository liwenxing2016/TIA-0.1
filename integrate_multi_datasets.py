#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import time
import sys

my_args = sys.argv

args_len = len(my_args)

fileCount = 0

exprMatrixPath = my_args[1]

datasetsInfo = []
for i in range(2, args_len):
    if my_args[i][0] == '[' and my_args[i][-1] == ']':
        datasetsInfo.append((my_args[i].split(',')[0][1:], my_args[i].split(',')[1][:-1]))
        fileCount = fileCount + 1

platformPath = my_args[fileCount+2]

samplePhenotypePath = my_args[fileCount+3]

resultsPath = my_args[fileCount+4]

if my_args[fileCount+5] == "FALSE":
    isAnnotated = False
else:
    isAnnotated = True

if my_args[fileCount+6] == "FALSE":
    isMerged = False
else:
    isMerged = True

if my_args[fileCount+7] == "FALSE":
    isIntegrate = False
else:
    isIntegrate = True

if my_args[fileCount+8] == "FALSE":
    isRenormalized = False
else:
    isRenormalized = True

merge_method = my_args[fileCount+9]

repeat_symbol = int(my_args[fileCount+10])

GSE_Files = []
platforms = []
for i in range(0, fileCount):
    GSE_Files.append(datasetsInfo[i][0])
    platforms.append(datasetsInfo[i][1])

annoExprMatrixPath = os.path.join(resultsPath, "2. Annotated Expression Matrix")
mergedExprMatrixPath = os.path.join(resultsPath, "3. Merged Expression Matrix")
integratedPath = os.path.join(resultsPath, "4. Integrated Dataset")
normalizedPath = os.path.join(resultsPath, "5. Normalized Dataset")

class Sequence:
    def __init__(self, sequence):
        # sequence of numbers we will process
        # convert all items to floats for numerical processing
        self.sequence = [float(item) for item in sequence]
    def count(self):
        return len(self.sequence)
    def sum(self):
        if self.count() < 1:
            return None
        else:
            return sum(self.sequence)
    def max(self):
        if self.count() < 1:
            return None
        else:
            return max(self.sequence)
    def min(self):
        if self.count() < 1:
            return None
        else:
            return min(self.sequence)
    def mean(self):
        if self.count() < 1:
            return None
        else:
            return sum(self.sequence) / self.count()
    def median(self):
        if self.count() < 1:
            return None
        else:
            self.sequence.sort()
            i = self.count() / 2
            if self.count() % 2 == 1:
                return self.sequence[i]
            else:
                return (self.sequence[i-1]+self.sequence[i]) / 2.0
    def std(self):
        if self.count() < 1:
            return None
        else:
            mean = self.mean()
            sq = sum([(i - mean) ** 2 for i in self.sequence])
            std = math.log((sq / (self.count() - 1)), 2)
            return std
    def gmean(self):
        if self.count() < 1:
            return None
        else:
            log_sum = 0
            for i in range(0, self.count()):
                log_sum = log_sum + math.log(self.sequence[i])
            return math.e ** (log_sum / self.count())
    def percentile(self, percentile):
        if self.count() < 1:
            value = None
        elif (percentile >= 100):
            sys.stderr.write("ERROR: percentile must be < 100. you supplied: %s\n" % percentile)
            value = None
        else:
            i = int(self.count() * (percentile / 100.0))
            self.sequence.sort()
            value = self.sequence[i]
        return value

def gene_annotation():
    print("Start gene annotation:")
    if not os.path.isdir(annoExprMatrixPath):
        os.mkdir(annoExprMatrixPath)
    for i in range(0, fileCount):
        f = open(os.path.join(exprMatrixPath, GSE_Files[i]), 'r')
        title = f.readline()
        exprMatrix = f.readlines()
        f.close()
        f = open(os.path.join(platformPath, platforms[i]), 'r')
        annoInfo = f.readlines()
        f.close()
        skip_line = 0
        probe_line = 0
        symbol_line = 0
        anno_type = 0
        for j in range(0, len(annoInfo)):
            if annoInfo[j][0] == '#':
                if "#ID" in annoInfo[j]:
                    probe_line = j
                if "#Gene Symbol" in annoInfo[j]:
                    symbol_line = j
                    anno_type = 0
                elif "#GeneSymbol" in annoInfo[j]:
                    symbol_line = j
                    anno_type = 0
                elif "#GENE_SYMBOL" in annoInfo[j]:
                    symbol_line = j
                    anno_type = 0
                elif "#Symbol" in annoInfo[j]:
                    symbol_line = j
                    anno_type = 0
                elif "#SYMBOL" in annoInfo[j]:
                    symbol_line = j
                    anno_type = 0
                elif "#ILMN_Gene" in annoInfo[j]:
                    symbol_line = j
                    anno_type = 0
                elif "#gene_assignment" in annoInfo[j]:
                    symbol_line = j
                    anno_type = 1
            else:
                skip_line = j + 1
                break
        annoInfo = annoInfo[skip_line:]
        annoInfo.sort()
        anno_len = len(annoInfo)
        annoExprMatrix = []
        annoExprMatrix.append(title)
        for line in exprMatrix:
            expr_probe = line.split('\t')[0]
            a = 0
            b = anno_len - 1
            count = 0
            while a <= b:
                j = int((a + b) / 2)
                anno_probe = annoInfo[j].split('\t')[probe_line]
                if expr_probe == anno_probe:
                    if anno_type == 0:
                        symbol = annoInfo[j].split('\t')[symbol_line]
                    elif anno_type == 1:
                        gene_assignment = annoInfo[j].split('\t')[symbol_line].split(" // ")
                        if len(gene_assignment) > 2:
                            symbol = gene_assignment[1]
                        else:
                            symbol = ""
                    if (symbol != "") and ("///" not in symbol) and ("---" not in symbol):
                        annoExprMatrix.append(symbol+line[len(expr_probe):])
                    break
                elif expr_probe > anno_probe:
                    a = j + 1
                elif expr_probe < anno_probe:
                    b = j - 1
                count = count + 1
                if count > math.log(anno_len, 2):
                    break
        f = open(os.path.join(annoExprMatrixPath, GSE_Files[i]), 'w')
        f.writelines(annoExprMatrix)
        f.close()
        print("Annotated expression matrix of %s has been writen to file." % GSE_Files[i])

def gene_symbol_merge(method="mean"):
    print("Start gene symbol merging, use \"%s\" method:" % method)
    if not os.path.isdir(mergedExprMatrixPath):
        os.mkdir(mergedExprMatrixPath)
    for i in range(0, fileCount):
        f = open(os.path.join(annoExprMatrixPath, GSE_Files[i]), 'r')
        title = f.readline()
        annoExprMatrix = f.readlines()
        f.close()
        annoExprMatrix.sort()
        matrix_row = len(annoExprMatrix)
        matrix_col = len(annoExprMatrix[0].split('\t'))
        mergedExprMatrix = []
        mergedExprMatrix.append(title)
        expr_i = 0
        while expr_i < matrix_row:
            symbol = annoExprMatrix[expr_i].split('\t')[0]
            tempMatrix = []
            for j in range(expr_i, matrix_row):
                if annoExprMatrix[j].split('\t')[0] == symbol:
                    tempMatrix.append(annoExprMatrix[j].strip().split('\t'))
                    expr_j = j
                else:
                    break
            tempMatrix = map(list, zip(*tempMatrix))
            tempLine = symbol
            for j in range(1, matrix_col):
                if method == "mean":
                    tempLine = tempLine + '\t' + str('%.6g' % Sequence(tempMatrix[j]).mean())
                elif method == "median":
                    tempLine = tempLine + '\t' + str('%.6g' % Sequence(tempMatrix[j]).median())
                elif method == "gmean":
                    tempLine = tempLine + '\t' + str('%.6g' % Sequence(tempMatrix[j]).gmean())
                elif method == "max":
                    tempLine = tempLine + '\t' + str('%.6g' % Sequence(tempMatrix[j]).max())
                elif method == "min":
                    tempLine = tempLine + '\t' + str('%.6g' % Sequence(tempMatrix[j]).min())
            mergedExprMatrix.append(tempLine + '\n')
            expr_i = expr_j + 1
        f = open(os.path.join(mergedExprMatrixPath, GSE_Files[i]), 'w')
        f.writelines(mergedExprMatrix)
        f.close()
        print("Merged expression matrix of %s has been writen to file." % GSE_Files[i])

def integrate_datasets(repeat_symbol=1):
    from collections import Counter
    temp_list = []
    print("Start to create total gene symbol list:")
    if repeat_symbol > fileCount:
        sys.stderr.write("ERROR: repeat symbol count must be < fileCount(%s). you supplied: %s\n" % (fileCount, repeat_symbol))
        return
    for i in range(0, fileCount):
        f = open(os.path.join(mergedExprMatrixPath, GSE_Files[i]), 'r')
        title = f.readline()
        mergedExprMatrix = f.readlines()
        f.close()
        for line in mergedExprMatrix:
            symbol = line.split('\t')[0]
            temp_list.append(symbol)
        print("File %s has been completed." % GSE_Files[i])
    symbol_dic = dict(Counter(temp_list))
    symbol_list = []
    for key in symbol_dic.keys():
        if symbol_dic[key] >= repeat_symbol:
            symbol_list.append(key)
    symbol_list.sort()
    print("The total number of genes: %d" % len(symbol_list))
    print("Start to integrate %d datasets:" % fileCount)
    if not os.path.isdir(integratedPath):
        os.mkdir(integratedPath)
    integrateTitle = "Symbol"
    integrateExprMatrix = []
    for symbol in symbol_list:
        integrateExprMatrix.append(symbol)
    for i in range(0, fileCount):
        f = open(os.path.join(mergedExprMatrixPath, GSE_Files[i]), 'r')
        integrateTitle = integrateTitle + '\t' + f.readline().strip()
        mergedExprMatrix = f.readlines()
        f.close()
        for line in mergedExprMatrix:
            symbol = line.split('\t')[0]
            a = 0
            b = len(symbol_list) - 1
            count = 0
            while a <= b:
                j = (a + b) / 2
                if symbol == symbol_list[j]:
                    integrateExprMatrix[j] = integrateExprMatrix[j] + line.strip()[len(symbol):]
                    break
                elif symbol > symbol_list[j]:
                    a = j + 1
                elif symbol < symbol_list[j]:
                    b = j - 1
                count = count + 1
                if count > math.log(len(symbol_list), 2):
                    break
        matrix_col = len(integrateTitle.split('\t'))
        for j in range(0, len(integrateExprMatrix)):
            current_col = len(integrateExprMatrix[j].split('\t'))
            for k in range(current_col, matrix_col):
                integrateExprMatrix[j] = integrateExprMatrix[j] + '\t' + "NA"
        print("File %s has been completed." % GSE_Files[i])
    integrateTitle = integrateTitle[7:] + '\n'
    for i in range(0, len(integrateExprMatrix)):
        integrateExprMatrix[i] = integrateExprMatrix[i] + '\n'
    f = open(os.path.join(integratedPath, "Integrated Dataset.txt"), 'w')
    f.write(integrateTitle)
    f.writelines(integrateExprMatrix)
    f.close()
    print("Integrated expression matrix has been writen to file.")

def global_renormalization():
    print("Start to global renormalization datasets:")
    if not os.path.isdir(normalizedPath):
        os.mkdir(normalizedPath)
    f = open(os.path.join(integratedPath, "Integrated Dataset.txt"), 'r')
    integrateTitle = f.readline()
    integrateExprMatrix = f.readlines()
    f.close()
    f = open(samplePhenotypePath, 'r')
    phenotypeTitle = f.readline()
    for i in range(0, len(phenotypeTitle.split('\t'))):
        if "sample" in phenotypeTitle.split('\t')[i].lower():
            sample_line = i
        if "study" in phenotypeTitle.split('\t')[i].lower():
            study_line = i
        if "category" in phenotypeTitle.split('\t')[i].lower():
            category_line = i
    select_study = []
    for i in range(0, len(datasetsInfo)):
        select_study.append(datasetsInfo[i][0])
    samplePhenotype = []
    for line in f.readlines():
        study = line.strip().split('\t')[study_line]
        if study in select_study:
            samplePhenotype.append(line)
    f.close()

    # map sample phenotype to integrated expression matrix
    # if sample phenotype haven't been sort, it can not use the Binary Search
    print("Map sample phenotype to integrated expression matrix.")
    tempSample = integrateTitle.strip().split('\t')
    tempPhenotype = []
    for i in range(0, len(samplePhenotype)):
        tempPhenotype.append(samplePhenotype[i])
    for i in range(0, len(tempSample)):
        for j in range(0, len(tempPhenotype)):
            if tempSample[i] == tempPhenotype[j].strip().split('\t')[sample_line]:
                samplePhenotype[i] = tempPhenotype[j]

    # locating sample, study IDs and category in the phenotype file
    sample = []
    study = []
    category = []
    for line in samplePhenotype:
        sample.append(line.strip().split('\t')[sample_line])
        study.append(line.strip().split('\t')[study_line])
        category.append(line.strip().split('\t')[category_line])

    dic_index = {}
    dic_cindex = {} # 'c' means control
    dic_means = {}
    dic_radio = {}
    value_list = []
    fkeys = GSE_Files

    # index of all sample in each study
    print("Calculate index of all sample in each study.")
    for i in range(0, len(study)):
        if study[i] not in dic_index.keys():
            dic_index[study[i]] = []
            dic_index[study[i]].append(i)
        else:
            dic_index[study[i]].append(i)

    # index of control sample in each study
    print("Calculate index of control sample in each study.")
    for i in range(0, len(study)):
        if category[i].lower() == "control":
            if study[i] not in dic_cindex.keys():
                dic_cindex[study[i]] = []
                dic_cindex[study[i]].append(i)
            else:
                dic_cindex[study[i]].append(i)

    # means of gene expression value in each study
    print("Calculate means of gene expression value in each study.")
    for line in integrateExprMatrix:
        symbol = line.split('\t')[0]
        line = line[len(symbol)+1:-1].split('\t')
        t_sum = 0 # 't' means total study
        t_len = 0 # 't' means total study
        for key in dic_cindex.keys():
            s_sum = 0 # 's' means a study
            s_len = 0 # 's' means a study
            for i in dic_cindex[key]:
                if line[i] != "NA" and line[i] != '':
                    isnull = False
                    s_sum = s_sum + 2 ** float(line[i])
                    s_len = s_len + 1
                    t_sum = t_sum + 2 ** float(line[i])
                    t_len = t_len + 1
                else:
                    isnull = True
            if key not in dic_means.keys():
                dic_means[key] = []
                if isnull == False:
                    dic_means[key].append('%.6g' % (s_sum / s_len))
                else:
                    dic_means[key].append('')
            else:
                if isnull == False:
                    dic_means[key].append('%.6g' % (s_sum / s_len))
                else:
                    dic_means[key].append('')
        if "All" not in dic_means.keys():
            dic_means["All"] = []
            dic_means["All"].append('%.6g' % (t_sum / t_len))
        else:
            dic_means["All"].append('%.6g' % (t_sum / t_len))

    # radio of gene expression value in each sample
    print("Calculate radio of gene expression value in each sample.")
    for key in dic_means.keys():
        for i in range(0, len(dic_means[key])):
            value = dic_means[key][i]
            ref_val = float(dic_means["All"][i])
            if key not in dic_radio.keys():
                dic_radio[key] = []
                if value != "NA" and value != '':
                    dic_radio[key].append('%.6g' % (ref_val / float(value)))
                else:
                    dic_radio[key].append('')
            else:
                if value != "NA" and value != '':
                    dic_radio[key].append('%.6g' % (ref_val / float(value)))
                else:
                    dic_radio[key].append('')

    # renormalization of each gene expression value
    print("Re-normalization of each gene expression value.")
    normalizeExprMatrix = []
    normalizeExprMatrix.append(integrateTitle)
    for i in range(0, len(integrateExprMatrix)):
        symbol = integrateExprMatrix[i].split('\t')[0]
        line = symbol
        value = integrateExprMatrix[i][len(symbol)+1:-1].split('\t')
        for key in fkeys:
            for j in dic_index[key]:
                if value[j] != "NA" and value[j] != '':
                    temp = '%.6g' % math.log((2 ** float(value[j])) * float(dic_radio[key][i]), 2)
                else:
                    temp = "NA"
                line = line + '\t' + str(temp)
        normalizeExprMatrix.append(line + '\n')

    f = open(os.path.join(normalizedPath, "Normalized Dataset.txt"), 'w')
    f.writelines(normalizeExprMatrix)
    f.close()
    print("Normalized expression matrix has been writen to file.")

if __name__ == "__main__":
    start = time.time()
    if isAnnotated == False:
        gene_annotation()
    if isMerged == False:
        gene_symbol_merge(merge_method)
    if isIntegrate == False:
        integrate_datasets(repeat_symbol)
    if isRenormalized == False:
        global_renormalization()
    end = time.time()
    print("Total Running Time: " + str(round(end-start, 3)) + "s.")
