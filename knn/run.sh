#!/bin/sh

DATA_DIR='data'

get_diabetes_data()
{
	wget -P $DATA_DIR https://archive.ics.uci.edu/ml/machine-learning-databases/00296/dataset_diabetes.zip	
	unzip data/dataset_diabetes.zip
	mv dataset_diabetes/* data/
	rm -rf dataset_diabetes
}

build_software()
{
	ant clean
	ant
}

run_preprocess()
{
	ant DiabetesMLPrep
	ant DiabetesMLNomn
}

benchmark()
{
	ant DiabetesMLBench
}

clusterize()
{
	ant DiabetesMLClust
}

reduce()
{
	ant DiabetesMLReduce
}

clusterize2()
{
	ant DiabetesMLClust2
}

############ P R O C E S S I N G ############
#get_diabetes_data
build_software
#run_preprocess
#benchmark
#clusterize
#reduce
clusterize2

##
