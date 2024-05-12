.. _inputFile:

##########
Input File
##########

The general parsing of the input file is based on a Lark grammar implementation (For more details see `Lark Grammar <https://lark-parser.readthedocs.io>`_). Any input file must be based on the following definitions of input key and value pairs:

.. note::
   There are two different types of input key and value pairs. The first type is the key and value pairs that are defined in line seperated by a :code:`=` e.g:

    .. code-block:: bash
    
        key = value

    The second type are so called multiline statements where in the first line the key is defined and in the following lines the values assigned to the key. The multiline statements must be closed by an :code:`END` statement. The following example shows a multiline statement:

    .. code-block:: bash

        key
        value1
        value2
        END

    It is important to note that multiline statements are always parsed as list/array like values. This means, if the documentation of the key states that the value is not a list or array, an inlined statement must be used.

.. note::
    In general, all keys are case-insensitive as well as the closing statement :code:`END` of a multiline statement. The values are case-sensitive. Furthermore, all keys and values are stripped from leading and trailing whitespaces and :code:`#` can be used to include comments (including inline comments). Inline statements using :code:`key = value` can also be used multiple times in one line separated by a :code:`,` to define multiple key and value pairs in one line e.g.:

    .. code-block:: bash

        key1 = value1, key2 = value2

.. note::
    The values are read as strings and are converted to the correct type based on the documentation of the key (if possible). In general, the user should not worry about the type of the value as the parser will try to convert the value to the correct type. If the conversion fails, an error will be raised. The following examples show the conversion of the values:

    * :code:`True` and :code:`False` are converted to :code:`bool` (case-insensitive) 
    * :code:`1` is converted to :code:`int` following possible conversions to :code:`float`
    * :code:`1.0` is converted to :code:`float` 
    * :code:`any-kind-of_string` is converted to :code:`str`
    * :code:`[1, 2, 3]` is converted to :code:`list` following possible (all values have to be of the   same type)
    * :code:`1..4` is converted to :code:`range` range(1, 4)
    * :code:`1-4` same as :code:`1..4` 
    * :code:`1..3..10` is converted to :code:`range` range(1, 10, 3), please note that the step size is always the middle value in contrast to the python syntax
    * :code:`1-3-10` same as :code:`1..3..10`
    * :code:`file_0*.text` is treated as a list of files matching the pattern (For more details see the `glob package <https://docs.python.org/3/library/glob.html>`_)




