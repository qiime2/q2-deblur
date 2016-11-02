import subprocess
import tempfile
import os
import qiime.plugin
import biom
from q2_types import FeatureTable, Frequency
from q2_types.per_sample_sequences import \
        SingleLanePerSampleSingleEndFastqDirFmt, FastqGzFormat
from q2_types.per_sample_sequences import SequencesWithQuality
from q2_types.feature_data import DNAIterator
from q2_types.sample_data import SampleData
from q2_types.feature_data import FeatureData, Sequence
import q2_deblur
from deblur.deblurring import get_default_error_profile


def workflow(demultiplexed_seqs: SingleLanePerSampleSingleEndFastqDirFmt,
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
             keep_tmp_files: bool=False,
             log_level: int=2,
             jobs_to_start: int=1) -> (biom.Table, DNAIterator):

    if error_dist is None:
        error_dist = get_default_error_profile()

    with tempfile.TemporaryDirectory() as tmp:
        iter_view = demultiplexed_seqs.sequences.iter_views(FastqGzFormat)

        for path, view in iter_view:
            print(path)
            #_single_sample(str(view), threads, tmp)


plugin = qiime.plugin.Plugin(
    name='q2-deblur',
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
    citation_text=None
)

# The next two code blocks are examples of how to register methods and
# visualizers. Replace them with your own registrations when you are ready to
# develop your plugin.

plugin.methods.register_function(
    function=workflow,
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
        'keep_tmp_files': qiime.plugin.Bool,
        'log_level': qiime.plugin.Int,
        'jobs_to_start': qiime.plugin.Int
    },
    outputs=[('table', FeatureTable[Frequency]),
             ('representative_sequences', FeatureData[Sequence])],
    name='Deblur',
    description='This method applies the Deblur workflow'
)
