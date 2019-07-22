import docker


def predict_protein_abundances(
        rna,
        cna,
        output_dir,
        tumor,
        logging=True,
        ):

    image_name = 'cptacdream/sub2:{}'.format(tumor)
    client = docker.from_env()

    if logging:
        print("Pulling image. This may take a few minutes...")

    client.images.pull(image_name)

    running_container = client.containers.run(
        image_name,
        detach=True,
        volumes={
            rna: {
                'bind': '/rna.txt',
                'mode': 'rw'
            },
            cna: {
                'bind': '/cna.txt',
                'mode': 'rw'
            },
            output_dir: {
                'bind': '/output',
                'mode': 'rw'
            }
        }
    )

    if logging:
        for line in running_container.logs(stream=True):
            print(line.strip())

    prediction_output_f = '{}/prediction.tsv'.format(output_dir)

    return prediction_output_f


def predict_phospho(
        rna,
        protein,
        output_dir,
        tumor,
        logging=True,
        ):

    image_name = 'cptacdream/sub3:{}'.format(tumor)
    client = docker.from_env()

    if logging:
        print("Pulling image. This may take a few minutes...")

    client.images.pull(image_name)

    running_container = client.containers.run(
        image_name,
        detach=True,
        volumes={
            rna: {
                'bind': '/rna.txt',
                'mode': 'rw'
            },
            protein: {
                'bind': '/proteome.txt',
                'mode': 'rw'
            },
            output_dir: {
                'bind': '/output',
                'mode': 'rw'
            }
        }
    )

    if logging:
        for line in running_container.logs(stream=True):
            print(line.strip())

    prediction_output_f = '{}/prediction.tsv'.format(output_dir)

    return prediction_output_f


if __name__ == '__main__':
    _container = predict_phospho(
        tumor='ovarian',
        rna='/Users/anna/Documents/DREAM_Challenge/hongyang_image_files/sub3_ova_CPTAC_submission/sub3/data/raw/retrospective_ova_rna_seq_sort_common_gene_15121.txt',
        protein='/Users/anna/Documents/DREAM_Challenge/hongyang_image_files/sub3_ova_CPTAC_submission/sub3/data/raw/retrospective_ova_PNNL_proteome_sort_common_gene_7061.txt',
        output_dir='/Users/anna/PycharmProjects/proteo_estimator/tests/output_breast_sub3',
        logging=True
        )
