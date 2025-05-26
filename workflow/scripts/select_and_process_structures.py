import os

import pandas as pd
from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.SeqUtils import IUPACData
from snakemake.script import snakemake
from stcrpy import fetch_TCRs

stcrdab_summary = pd.read_csv(snakemake.input[0], delimiter='\t')
stcrdab_summary['resolution'] = pd.to_numeric(stcrdab_summary['resolution'], errors='coerce')

print('Number of structures', len(stcrdab_summary))

selected_stcrdab = stcrdab_summary.query(
    "TCRtype == 'abTCR' and resolution <= 3.5 and antigen_type == 'peptide' and mhc_type in ('MH1', 'MH2')"
).copy()

print('Number of structures', len(selected_stcrdab))

output_dir = snakemake.output[0]

if not os.path.exists(output_dir):
    os.mkdir(output_dir)


selected_structures_summary = []

for pdb_id in selected_stcrdab['pdb'].unique():
    print('Loading', pdb_id)

    try:
        tcrs = fetch_TCRs(pdb_id)

    except PDBConstructionException as err:
        print('Skipping', pdb_id, 'because of', err)
        continue

    for tcr in tcrs:
        tcr_info = {
            'pdb_id': pdb_id,
            'alpha_chain_id': tcr.VA,
            'beta_chain_id': tcr.VB,
        }

        for cdr in tcr.get_CDRs():
            tcr_info[cdr.id + '_sequence'] = ''.join([
                IUPACData.protein_letters_3to1[res.get_resname().title()] for res in cdr
                if res.get_resname().title() in IUPACData.protein_letters_3to1
            ])

        if len(tcr.antigen) > 1:
            print('Discarding multi-antigen')
            continue

        try:
            tcr_info['antigen_chain_id'] = tcr.antigen[0].id
        except IndexError as err:
                print('Skipping', pdb_id, 'becasue of', err)
                continue

        tcr_info['antigen_sequence'] = ''.join([
            IUPACData.protein_letters_3to1[res.get_resname().title()] for res in tcr.antigen[0]
            if res.get_resname().title() in IUPACData.protein_letters_3to1
        ])

        try:
            entry_sequences = (
                tcr_info['cdra1_sequence'],
                tcr_info['cdra2_sequence'],
                tcr_info['cdra3_sequence'],
                tcr_info['cdrb1_sequence'],
                tcr_info['cdrb2_sequence'],
                tcr_info['cdrb3_sequence'],
                tcr_info['antigen_sequence'],
            )
        except KeyError as err:
            print('Skipping', pdb_id, 'becasue of', err)
            continue

        already_added = False
        for past_entry in selected_structures_summary:
            past_entry_sequences = (
                past_entry['cdra1_sequence'],
                past_entry['cdra2_sequence'],
                past_entry['cdra3_sequence'],
                past_entry['cdrb1_sequence'],
                past_entry['cdrb2_sequence'],
                past_entry['cdrb3_sequence'],
                past_entry['antigen_sequence'],
            )

            if entry_sequences == past_entry_sequences:
                already_added = True
                break

        if already_added:
            print('Skipping TCR with same CDR and antigen sequences')
            continue

        if len(tcr.MHC) > 1:
            print('Discarding multi-MHC')
            continue

        try:
            tcr_info['mhc1_chain_id'] = tcr.MHC[0].id[0]

        except IndexError as err:
            print('Skipping', pdb_id, 'because of', err)
            continue

        tcr_info['mhc2_chain_id'] = tcr.MHC[0].id[1] if len(tcr.MHC[0]) > 1 else None

        for key, info in tcr.get_germlines_and_alleles().items():
            match key:
                case 'VA':
                    tcr_info['v_alpha_gene'] = info[0]
                    tcr_info['j_alpha_gene'] = info[1]

                case 'VB':
                    tcr_info['v_beta_gene'] = info[0]
                    tcr_info['j_beta_gene'] = info[1]

                case 'VA_species':
                    tcr_info['alpha_chain_species'] = info

                case 'VB_species':
                    tcr_info['beta_chain_species'] = info

                case 'TCR_VA_seq':
                    tcr_info['alpha_chain_sequence'] = info

                case 'TCR_VB_seq':
                    tcr_info['beta_chain_sequence'] = info

                case 'MHC_GA1' | 'MHC_GA':
                    tcr_info['mhc1_gene'] = info

                case 'MHC_GA1_seq' | 'MHC_GA_seq':
                    tcr_info['mhc1_sequence'] = info

                case 'MHC_B2M' | 'MHC_GB':
                    tcr_info['mhc2_gene'] = info

                case 'MHC_B2M_seq' | 'MHC_GB_seq':
                    tcr_info['mhc2_sequence'] = info

                case 'antigen':
                    pass

                case _:
                    tcr_info[key.lower()] = info

        try:
            docking_geometry = tcr.calculate_docking_geometry(mode='cys')

        except (ValueError, ZeroDivisionError, KeyError) as err:
            print('Skipping', pdb_id, 'because of', err)
            continue

        tcr_info['pitch_angle'] = docking_geometry['pitch_angle']
        tcr_info['scanning_angle'] = docking_geometry['scanning_angle']

        tcr_name = (
            f"{pdb_id}_"
            f"{tcr_info['alpha_chain_id']}{tcr_info['beta_chain_id']}"
            f"{tcr_info['antigen_chain_id']}"
            f"{tcr_info['mhc1_chain_id']}{tcr_info['mhc2_chain_id']}"
        )

        tcr.standardise_chain_names()

        try:
            interactions = tcr.profile_peptide_interactions()

        except KeyError as err:
            print('Skipping structure', pdb_id, 'because of', err)
            continue

        interactions.to_csv(os.path.join(output_dir, tcr_name + '_interactions.csv'), index=False)
        tcr.save(os.path.join(output_dir, tcr_name + '.pdb'))

        selected_structures_summary.append(tcr_info)

selected_structures_summary = pd.DataFrame(selected_structures_summary)
selected_structures_summary

selected_structures_summary.to_csv(os.path.join(output_dir, 'structures_summary.csv'), index=False)
