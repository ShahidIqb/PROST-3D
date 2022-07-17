from six.moves import urllib

#This function is taken from https://github.com/oddt/oddt/tree/master/oddt
def fetchstructure(pdbid):
    """Fetch pdb structure using pdbid."""
    req = urllib.request.Request('https://files.rcsb.org/view/%s.pdb' % pdbid)
    response = urllib.request.urlopen(req)
    pdb_block = response.read().decode('utf-8')
    return pdb_block
