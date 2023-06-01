from upath import UPath
from upath.implementations.cloud import GCSPath

groups = [
    "20230317_06h00m45s",
    "20230410_21h06m43s",
]

path_root = GCSPath("gs://liulab/differential_composition_and_expression")

for group in groups:
    things = path_root.glob(f"{group}/**}")
    
    if f.is_dir():
        print(f)

path_root_2 = UPath("gs://liulab/differential_composition_and_expression/20230410_21h06m43s")
