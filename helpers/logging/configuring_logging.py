import logging

LOGGING_FORMAT = "%(asctime)s %(process)d/%(threadName)s %(name)s %(levelname)s %(message)s"
formatter = logging.Formatter(LOGGING_FORMAT)


def configure_logging():
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logging.getLogger().handlers = [handler]
