import subprocess
import tempfile
import os
import qiime.plugin
import biom
import skbio
import subprocess
import hashlib

from q2_types import FeatureTable, Frequency
from q2_types.per_sample_sequences import \
        SingleLanePerSampleSingleEndFastqDirFmt, FastqGzFormat
from q2_types.per_sample_sequences import SequencesWithQuality
from q2_types.feature_data import DNAIterator
from q2_types.sample_data import SampleData
from q2_types.feature_data import FeatureData, Sequence
import q2_deblur
from deblur.deblurring import get_default_error_profile


def denoise(demultiplexed_seqs: SingleLanePerSampleSingleEndFastqDirFmt,
            pos_ref_fp: str=None,
            neg_ref_fp: str=None,
            mean_error: float=0.005, 
            error_dist: str=None, 
            indel_prob: float=0.01, 
            indel_max: int=3,
            trim_length: int=100, 
            min_reads: int=0, 
            min_size: int=2, 
            negate: bool=False, 
            jobs_to_start: int=1,
            hashed_feature_ids: bool=True) -> (biom.Table, DNAIterator):

    if error_dist is None:
        error_dist = get_default_error_profile()

    with tempfile.TemporaryDirectory() as tmp:
        seqs_fp = str(demultiplexed_seqs)
        cmd = ['deblur', 'workflow', 
               '--seqs-fp', seqs_fp,
               '--output-dir', tmp,
               '--mean-error', str(mean_error),
               '--error-dist', ','.join([str(i) for i in error_dist]),
               '--indel-prob', str(indel_prob),
               '--indel-max', str(indel_max),
               '--trim-length', str(trim_length),
               '--min-reads', str(min_reads),
               '--min-size', str(min_size),
               '-w']  
        if pos_ref_fp is not None:
            cmd.append('--pos-ref-db')
            cmd.append(pos_ref_fp)
        
        if neg_ref_fp is not None:
            cmd.append('--neg-ref-db')
            cmd.append(neg_ref_fp)
        
        if negate:
            cmd.append('--negate')
        
        subprocess.run(cmd, check=True)
        
        # code adapted from q2-dada2
        table = biom.load_table(os.path.join(tmp, 'final.biom'))
        sid_map = {id_: id_.split('_')[0] for id_ in table.ids(axis='sample')}
        table.update_ids(sid_map, axis='sample', inplace=True)
        
        if hashed_feature_ids:
            # Make feature IDs the md5 sums of the sequences.
            fid_map = {id_: hashlib.md5(id_.encode('utf-8')).hexdigest()
                       for id_ in table.ids(axis='observation')}
            table.update_ids(fid_map, axis='observation', inplace=True)

            rep_sequences = DNAIterator((skbio.DNA(k, metadata={'id': v},
                                                   lowercase='ignore')
                                         for k, v in fid_map.items()))
        else:
            rep_sequences = DNAIterator(
                (skbio.DNA(id_, metadata={'id': id_}, lowercase='ignore')
                 for id_ in table.ids(axis='observation')))

    return (table, rep_sequences)


plugin = qiime.plugin.Plugin(
    name='deblur',
    version=q2_deblur.__version__,
    website='https://github.com/biocore/deblur',
    package='q2_deblur',
    # Information on how to obtain user support should be provided as a free
    # text string via user_support_text. If None is provided, users will
    # be referred to the plugin's website for support.
    user_support_text='https://github.com/biocore/deblur/issues',
    # Information on how the plugin should be cited should be provided as a
    # free text string via citation_text. If None is provided, users
    # will be told to use the plugin's website as a citation.
    citation_text=("Deblur rapidly resolves single-nucleotide community "
                   "sequence patterns. Amnon Amir, Daniel McDonald, Jose "
                   "A. Navas-Molina, Evguenia Kopylova, Jamie Morton, "
                   "Zhenjiang Zech Xu, Eric P. Kightley, Luke R. Thompson, "
                   "Embriette R. Hyde, Antonio Gonzalez, Rob Knight. mSystems "
                   "(in press).")
)

plugin.methods.register_function(
    function=denoise,
    inputs={
        'demultiplexed_seqs': SampleData[SequencesWithQuality]
    },
    parameters={
        'pos_ref_fp': qiime.plugin.Str,
        'neg_ref_fp': qiime.plugin.Str,
        'mean_error': qiime.plugin.Float,
        'error_dist': qiime.plugin.Str,
        'indel_prob': qiime.plugin.Float,
        'indel_max': qiime.plugin.Int,
        'trim_length': qiime.plugin.Int,
        'min_reads': qiime.plugin.Int,
        'min_size': qiime.plugin.Int,
        'negate': qiime.plugin.Bool,
        'jobs_to_start': qiime.plugin.Int,
        'hashed_feature_ids': qiime.plugin.Bool
    },
    outputs=[('table', FeatureTable[Frequency]),
             ('representative_sequences', FeatureData[Sequence])],
    name='Deblur',
    description='This method applies the Deblur workflow'
)
