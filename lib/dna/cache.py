# cache.py - Caches values from Genbank in buckets so as not to duplicate calls
# Stores gene IDs, gene sequences, etc., separately

import os
import yaml

cache = { 
    'file': 'Cache.yml',
    'genes': {},          # name -> gene ID (in 'gene')
    'coordinates': {},    # gene ID -> (gi, strand, start, end)
    'sequences': {}       # (gi, strand, start, end) -> sequence
}

def load_cache(root=os.getcwd()):
    global cache
    cache_file = root + cache['file']

    if os.path.isfile(cache_file):
        with open(cache_file, 'r') as f:
            cache = yaml.load(f)

def save_cache(root=os.getcwd()):
    global cache
    cache_file = root + cache['file']
    with open(cache_file, 'w') as f:
        f.write(yaml.dump(cache))
