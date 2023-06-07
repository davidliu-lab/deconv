import os
from pathlib import Path

from dagster import get_dagster_logger, job, op


@op
def get_file_sizes():
    # files = [f for f in os.listdir("/home/jupyter/deconv") if os.path.isfile(f)]
    files = [p for p in Path("/home/jupyter/deconv").iterdir() if p.is_file()]
    return {p: p.stat().st_size for p in files}


@op
def report_total_size(file_sizes):
    total_size = sum(file_sizes.values())
    # In real life, we'd send an email or Slack message instead of just logging:
    get_dagster_logger().info(f"Total size: {total_size}")


@job
def serial():
    report_total_size(get_file_sizes())


@op
def get_total_size(file_sizes):
    return sum(file_sizes.values())


def test_get_total_size():
    file_sizes = {"foo": 1, "bar": 2}
    result = get_total_size(file_sizes)
    assert result == 3


@op
def get_largest_size(file_sizes):
    return max(file_sizes.values())


@op
def report_file_stats(total_size, largest_size):
    # In real life, we'd send an email or Slack message instead of just logging:
    get_dagster_logger().info(f"Total size: {total_size}, largest size: {largest_size}")


@job
def diamond():
    file_sizes = get_file_sizes()
    total_size = get_total_size(file_sizes)
    print(total_size)
    report_file_stats(
        total_size=total_size,
        largest_size=get_largest_size(file_sizes),
    )


def test_diamond():
    result = diamond.execute_in_process()
    assert result.success
    assert result.output_for_node("get_total_size") > 0


if __name__ == "__main__":
    result = diamond.execute_in_process()
    print(result)
