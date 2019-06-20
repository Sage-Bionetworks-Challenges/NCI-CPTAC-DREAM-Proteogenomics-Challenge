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


if __name__ == '__main__':
    _container = predict_protein_abundances(
        tumor='ovarian',
        rna='/Users/anna/Documents/DREAM_Challenge/hongyang_image_files/sub2_breast_CPTAC_breast/rna.txt',
        cna='/Users/anna/Documents/DREAM_Challenge/hongyang_image_files/sub2_breast_CPTAC_breast/cna.txt',
        output_dir='/Users/anna/PycharmProjects/proteo_estimator/tests/output_ova',
        logging=True
        )
