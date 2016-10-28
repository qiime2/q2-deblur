import qiime.plugin

import q2_deblur

# These imports are only included to support the example methods and
# visualizers. Remove these imports when you are ready to develop your plugin.
from q2_dummy_types import IntSequence1, IntSequence2, Mapping
from ._dummy_method import concatenate_ints
from ._dummy_visualizer import mapping_viz

plugin = qiime.plugin.Plugin(
    name='q2-deblur',
    version=q2_deblur.__version__,
    website='https://github.com/biocore/deblur',
    package='q2_deblur',
    # Information on how to obtain user support should be provided as a free
    # text string via user_support_text. If None is provided, users will
    # be referred to the plugin's website for support.
    user_support_text=None,
    # Information on how the plugin should be cited should be provided as a
    # free text string via citation_text. If None is provided, users
    # will be told to use the plugin's website as a citation.
    citation_text=None
)

# The next two code blocks are examples of how to register methods and
# visualizers. Replace them with your own registrations when you are ready to
# develop your plugin.

plugin.methods.register_function(
    function=concatenate_ints,
    inputs={
        'ints1': IntSequence1 | IntSequence2,
        'ints2': IntSequence1,
        'ints3': IntSequence2
    },
    parameters={
        'int1': qiime.plugin.Int,
        'int2': qiime.plugin.Int
    },
    outputs=[
        ('concatenated_ints', IntSequence1)
    ],
    name='Concatenate integers',
    description='This method concatenates integers into a single sequence in '
                'the order they are provided.'
)

plugin.visualizers.register_function(
    function=mapping_viz,
    inputs={
        'mapping1': Mapping,
        'mapping2': Mapping
    },
    parameters={
        'key_label': qiime.plugin.Str,
        'value_label': qiime.plugin.Str
    },
    name='Visualize two mappings',
    description='This visualizer produces an HTML visualization of two '
                'key-value mappings, each sorted in alphabetical order by key.'
)
