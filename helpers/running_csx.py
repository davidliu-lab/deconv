import os


class InputFile:
    def __init__(self, target_filename, source_uri):
        self.target_filename = target_filename
        self.source_uri = source_uri

    @property
    def source_local_path(self):
        return self.source_uri.replace("gs://", "/mnt/buckets/")


class DockerJob:
    def __init__(self, path, csx_input_files, other_args):
        self.path = path
        self.csx_input_files = csx_input_files
        self.other_args = other_args

    def make_copy_commands(self):
        commands = list()
        for file_arg, input_file in self.csx_input_files.items():
            target_path = os.path.join(self.path, "in", input_file.target_filename)
            commands.append(f"gsutil cp {input_file.source_uri} {target_path}")
        return commands

    def make_docker_command(self):
        command = f"""docker run \
    --rm \
    -v {self.path}/in:/src/data \
    -v {self.path}:/src/outdir \
    --user "$(id -u):$(id -g)" \
    cibersortx/fractions:latest \
    --username lyronctk@stanford.edu \
    --token dfeba2c8b9d61daebee5fa87026b8e56 \
    --replicates 5 \
    --sampling 0.5 \
    --fraction 0.75 \
    --k.max 999 \
    --q.value 0.01 \
    --G.min 300 \
    --G.max 500 \
    --filter FALSE \
    --QN FALSE \\""".replace(
            "     ", " \\\n    "
        )
        for arg, input_file in self.csx_input_files.items():
            command += f"\n    --{arg} {input_file.target_filename} \\"
        for arg, value in self.other_args.items():
            command += f"\n    --{arg} {value} \\"
        command += "\n    --verbose TRUE"
        return command


class Experiment:
    def __init__(self, root_uri, name):
        self.root_uri = root_uri
        self.name = name

    @property
    def gcs_uri(self):
        return os.path.join(self.root_uri, self.name)

    @property
    def local_path(self):
        return self.gcs_uri.replace("gs://", "/mnt/buckets/")
