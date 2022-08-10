import datetime


def make_a_nice_timestamp_of_now() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
