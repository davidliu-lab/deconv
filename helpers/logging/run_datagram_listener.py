"""Script to run a datagram server, listening for logging messages and emitting them to the console"""

import asyncio
import logging
import pickle

from .configuring_logging import formatter


class DatagramLoggingProtocol(asyncio.DatagramProtocol):
    """Protocol to listen for logging messages and emit them to the console"""

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    def datagram_received(self, data, addr):
        obj = pickle.loads(data[4:])
        log_record = logging.makeLogRecord(obj)
        self.handler.emit(log_record)


if __name__ == "__main__":
    loop = asyncio.get_event_loop()
    listen = loop.create_datagram_endpoint(
        DatagramLoggingProtocol, local_addr=("127.0.0.1", 12000)
    )
    transport, protocol = loop.run_until_complete(listen)
    try:
        loop.run_forever()
    except KeyboardInterrupt:
        pass
    finally:
        print("Closing datagram transport")
        transport.close()
        loop.close()
