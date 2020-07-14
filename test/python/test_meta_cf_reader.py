from teca import *
import sys

set_stack_trace_on_error()

def parse_args(args):
    regex = []
    var = []
    out_file = ''
    baseline = ''
    time_reader = ''
    geometry_reader = ''
    tmp_vars = []
    in_list = False
    it = iter(args)
    while True:
        try:
            arg = next(it)

            if arg == ')':
                in_list = False
                var.append(tmp_vars)

            elif in_list:
                tmp_vars.append(arg)

            elif arg == '(':
                in_list = True
                tmp_vars = []
                try:
                    regex.append(next(it))
                except Exception:
                    raise RuntimeError('Missing regex at begining of list')

            elif arg == '-t':
                try:
                    time_reader = next(it)
                except Exception:
                     raise RuntimeError('Missing time reader key')

            elif arg == '-g':
                try:
                    geometry_reader = next(it)
                except Exception:
                     raise RuntimeError('Missing geometry reader key')

            elif arg == '-o':
                try:
                    out_file = next(it)
                except Exception:
                    raise RuntimeError('Missing out file name')

            elif arg == '-b':
                try:
                    baseline = next(it)
                except Exception:
                    raise RuntimeError('Missing base file name')

            elif arg == '-h':
                sys.stderr.write('usage: test_meta_cf_reader [-o out file] '
                                 '[-b base line] [( regex var0 ... varn )] '
                                 ' ... [( regex var0 ... varn )]\n')
                sys.exit(-1)

        except StopIteration:
            return regex, time_reader, geometry_reader, var, out_file, baseline


regex, time_reader, geometry_reader, var, out_file, baseline = parse_args(sys.argv)


print(regex)
print(time_reader)
print(geometry_reader)
print(var)
print(out_file)
print(baseline)



cfr = teca_meta_cf_reader.New()
n = len(regex)
i = 0
while i < n:
    key = 'r_%d'%(i)
    cfr.add_reader(key, regex[i], 0, 0, var[i])
    i += 1


cfr.set_x_axis_variable('lon')
cfr.set_y_axis_variable('lat')
cfr.set_t_axis_variable('time')
cfr.set_time_reader(time_reader)
cfr.set_geometry_reader(geometry_reader)

md = cfr.update_metadata()
sys.stderr.write('md = %s\n'%(str(md)))
sys.exit(0)




coords = teca_normalize_coordinates.New()
coords.set_input_connection(cfr.get_output_port())

exe = teca_index_executive.New()
exe.set_start_index(first_step)
exe.set_end_index(end_index)
exe.set_arrays(arrays)

wri = teca_cartesian_mesh_writer.New()
wri.set_input_connection(coords.get_output_port())
wri.set_executive(exe)
wri.set_file_name(out_file)

wri.update()

