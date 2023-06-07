from kubernetes import client, config
from kubernetes.client.rest import ApiException


def run_image(image, command, cpu, memory):
    # Load the Kubernetes configuration
    config.load_kube_config()

    # Create an instance of the Kubernetes API client
    api_instance = client.BatchV1Api()

    # Define the container specification
    container = client.V1Container(
        name="my-container",
        image=image,
        command=command.split(),
        resources=client.V1ResourceRequirements(
            limits={"cpu": cpu, "memory": memory}, requests={"cpu": cpu, "memory": memory}
        ),
    )

    # Define the pod specification
    pod_spec = client.V1PodSpec(containers=[container], restart_policy="Never")

    # Define the job specification
    job_spec = client.V1JobSpec(template=client.V1PodTemplateSpec(spec=pod_spec), backoff_limit=0)

    # Define the job
    job = client.V1Job(
        api_version="batch/v1",
        kind="Job",
        metadata=client.V1ObjectMeta(name="my-job"),
        spec=job_spec,
    )

    try:
        # Create the job
        api_instance.create_namespaced_job(namespace="default", body=job)

        # Wait for the job to complete
        api_instance.read_namespaced_job_status(namespace="default", name="my-job", watch=True)

        # Get the logs from the completed pod
        logs = api_instance.read_namespaced_pod_log(
            name="my-job-xxxxx", namespace="default"  # Replace with the actual pod name
        )

        # Process the logs as desired
        print(logs)

    except ApiException as e:
        print(f"Exception when creating or retrieving the job: {e}")


# Example usage
image = "grisaitis/bayesprism:latest"
command = "Rscript run_demo_gbm_simple.R"  # Replace x with your desired value
cpu = "5"
memory = "20Gi"

run_image(image, command, cpu, memory)
