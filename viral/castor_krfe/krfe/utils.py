from sklearn.externals import joblib

def save_model(clf, model_file):
    # Save model
    joblib.dump(clf, model_file) 


def load_model(model_file):
    # Load model
    clf = joblib.load(model_file)
 
    return clf


def fetch_list_from_file(my_file):
    with open(my_file) as fh:
        lines = fh.read().splitlines()

        return lines
