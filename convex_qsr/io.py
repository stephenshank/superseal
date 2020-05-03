import json

from Bio import SeqIO


def write_json(path, data):
    if path:
        with open(path, 'w') as json_file:
            json.dump(data, json_file, indent=2)


def write_fasta(path, data):
    if path:
        SeqIO.write(data, path, 'fasta')


def write_csv(path, data):
    if path:
        data.to_csv(path)
