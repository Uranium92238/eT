def configure(options, input_files, extra_args):
    """
    This function is used by runtest to configure runtest
    at runtime for code specific launch command and file naming.
    """

    from os import path

    launcher = 'eT_launch'
    launcher_full_path = path.normpath(path.join(options.binary_dir, launcher))

    inp = input_files

    command = []
    command.append('{0}'.format(launcher_full_path))
    command.append('{0} '.format(inp))
    if extra_args is not None:
        command.append(extra_args)

    full_command = ' '.join(command)

    inp_no_suffix = path.splitext(inp)[0]

    output_prefix = '{0}'.format(inp_no_suffix)

    relative_reference_path = 'result'

    return launcher, full_command, output_prefix, relative_reference_path
