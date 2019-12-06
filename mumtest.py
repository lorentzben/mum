import argparse
import datetime

# Build Argument Parser in order to facilitate ease of use for user
parser = argparse.ArgumentParser(
    description="Perform Analysis of 16s Microbial Data using Mothur")
parser.add_argument('-n', action='store', required=True,
                    help="name for analysis, will be filename for resultant files", dest='job_name')
parser.add_argument('-c', action='store_true', default=False, required=False,
                    help="flag to use custom parameters, ignore for kelly lab parameters", dest='custom')
parser.add_argument('-l', action='store', default=300, type=int, required=False,
                    help='determines the longest value permitted for sequences', dest='ml')
parser.add_argument('-p', action='store', default=2, type=int, required=False,
                    help='pre cluster value, higher is more stringent', dest='pre')
parser.add_argument('-s', action='store', type=int, required=False,
                    help='sub sample value, only applicable for custom runs', dest='sub')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')


args = parser.parse_args()
print(args)
# main()

#esults = parser.parse_args()
#print 'simple_value     =', results.simple_value
#print 'constant_value   =', results.constant_value
#print 'boolean_switch   =', results.boolean_switch
#print 'collection       =', results.collection
#print 'const_collection =', results.const_collection