"Simple input/output helper functions for various formats."
import json

from Bio import SeqIO


def write_json(path, data):
    "Write a dictionary to a JSON file."
    if path:
        with open(path, 'w') as json_file:
            json.dump(data, json_file, indent=2)


def write_fasta(path, data):
    "Write an array of sequence records to a FASTA file."
    if path:
        SeqIO.write(data, path, 'fasta')


def write_csv(path, data):
    "Write a pandas dataframe to a CSV file."
    if path:
        data.to_csv(path)


def read_json(path):
    "Read a JSON file into a dictionary."
    with open(path) as json_file:
        data = json.load(json_file)
    return data
