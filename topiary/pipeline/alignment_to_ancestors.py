
import topiary

def generate_ancestors(df,output_dir,do_bootstrap=True):
    """
    """

    topiary.find_best_model(df,
                            output="00_find-model")

    topiary.generate_ml_tree(previous_dir="00_find-model",
                             output="01_ml-tree",
                             bootstrap=do_bootstrap)

    topiary.reconcile(previous_dir="01_ml-tree",
                      output="02_reconciliation")

    topiary.generate_ancestors(previous_dir="02_reconciliation",
                               output="03_ancestors")
