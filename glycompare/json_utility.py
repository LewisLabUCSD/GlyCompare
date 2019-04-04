import json


def store_json(address, dic):
    with open(address, 'w') as fp:
        json.dump(dic, fp)


def load_json(address):
    with open(address, 'r') as f:
        return json.load(f)

