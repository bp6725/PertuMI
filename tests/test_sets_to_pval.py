from pipelines.setes_to_pval import SetsToPval
from pipelines.pipelines_utils import PipelinesUtils
from pipelines.pipelines_codebase import PipelinesCodebase
from DAL.builder_synthetic_data import BuilderSyntheticData

run_sets = PipelinesCodebase.build_synthetic_experiment_sets(3,0.1)
all_res = SetsToPval.run_pipeline(*run_sets, n_permut=100)
res = all_res['original'][0]
print("end")