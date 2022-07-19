import logging


format_string = (
    "%(asctime)s %(process)d/%(threadName)s %(name)s %(levelname)s\n%(message)s"
)
formatter = logging.Formatter(format_string)

def configure_logging():
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logging.getLogger().handlers = [handler]
