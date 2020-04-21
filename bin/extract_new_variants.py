import baltic as bt
import argparse
import os, json
import pandas as pd

def label_reference_subtree(ll, new_sample_string):
    """Labels the subtree of reference sequences, ie, those which are
    ancestors of non-new sequences."""

    for node in ll.Objects:
        node.traits['ref'] = False

    ref_nodes = ll.getExternal(lambda k: new_sample_string not in k.traits['name'])
    while len(ref_nodes) > 0:
        node = ref_nodes.pop()
        node.traits['ref'] = True
        if node.parent and not node.parent.traits.get('ref'):
            ref_nodes.append(node.parent)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='')
    parser.add_argument("--pipeline_dir")
    parser.add_argument("--out_path")
    parser.add_argument("--new_sample_string", default='RR057')
    args = parser.parse_args()

    with open(os.path.join(args.pipeline_dir, 'ncov.json')) as fp:
        nx_tree=json.load(fp)

    tree = nx_tree['tree']
    meta = nx_tree['meta']

    json_translation = {'absoluteTime': lambda k: k.traits['node_attrs']['num_date']['value'],'name':'name'}
    ll = bt.loadJSON(tree, json_translation)

    label_reference_subtree(ll, args.new_sample_string)
    new_leaves = ll.getExternal(lambda k: not k.traits['ref'])
    leaves = ll.getExternal()

    new_mutations = set()
    introductions = set()
    for node in new_leaves:
        node_muts = []
        n = node
        while not n.traits['ref']:
            if 'branch_attrs' in n.traits and 'nuc' in n.traits['branch_attrs']['mutations']:
                muts = n.traits['branch_attrs']['mutations']['nuc']
                for m in muts:
                    node_muts.append(m)
                    new_mutations.add(m)
            n = n.parent
        introductions.add(n)
        node.traits['node_attrs']['new_mutations'] = node_muts

    print(f"Number of clusters: {len(introductions)}")

    for node in leaves:
        node_muts = []
        n = node
        while n is not None:
            if 'branch_attrs' in n.traits and 'nuc' in n.traits['branch_attrs']['mutations']:
                muts = n.traits['branch_attrs']['mutations']['nuc']
                for m in muts:
                    #if m in new_mutations:
                    node_muts.append(m)
            n = n.parent
        node.traits['node_attrs']['mutations'] = node_muts

    mutations_df = []
    for node in leaves:
        name = node.traits['name']
        for mut in node.traits['node_attrs']['mutations']:
            try:
                new_mut = mut in node.traits['node_attrs']['new_mutations']
            except KeyError:
                new_mut = False
            info = {
                'sample': name,
                'pos': int(mut[1:-1]),
                'ref': mut[0],
                'allele': mut[-1],
                'type': 'SNP',
                'new_mutation': new_mut}
            mutations_df.append(info)
    df = pd.DataFrame(mutations_df)
    df = df[['sample', 'pos', 'ref', 'allele', 'type', 'new_mutation']].sort_values(['pos', 'sample'])
    df.to_csv(args.out_path, index=False)
