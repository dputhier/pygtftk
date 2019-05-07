import re
import textwrap
import xml.etree.ElementTree as ElementTree
from collections import defaultdict

import requests
from requests.exceptions import ConnectionError as ConErr

from pygtftk.utils import message


class Biomart(object):

    def __init__(self,
                 xml=None,
                 url="http://www.ensembl.org/biomart/martservice?",
                 http_proxy='',
                 https_proxy=''):
        """A connection to a Mart server.

        :param url: The full url to the mart server.
        :param xml: A string/docstring containing the query.
        :param http_proxy:

        """

        self.url = url
        self.xml = xml
        self.proxies = {'http': http_proxy, 'https': https_proxy}
        self.databases = []
        self.datasets = defaultdict(list)
        self._get_databases()
        self.response = None

    def __repr__(self, *args, **kwargs):

        msg = textwrap.dedent("""
        - A Biomart connection to {a} databases:

        """.format(a=len(self.databases)))

        for i in self.databases:
            msg += "    -" + str(i) + "\n"

        msg += textwrap.dedent("""
        - {n} datasets are available.

        """.format(n=len(self.datasets)))

        if len(self.datasets) == 0:
            msg += textwrap.dedent("""
                - Use the .get_data_sets() methods with one database as
                  argument (e.g .get_data_sets('{d}').)
            """.format(d=self.databases[0]))

        return msg

    def _get_databases(self):

        message("Listing available databases", type="DEBUG")
        try:
            self.query(query={'type': 'registry'})
        except ConErr:
            message("Raised a connection Error.", type="ERROR")

        tree = ElementTree.fromstring(self.response.content)

        for child in tree:
            if child.tag == 'MartURLLocation':
                self.databases += [child.attrib['name']]

    def get_datasets(self, database=None):
        message("Listing available datasets", type="DEBUG")
        if database in self.databases:

            self.query(query={'type': 'datasets', 'mart': database})
            for i in self.response.text.split("\n"):
                fields = i.split("\t")
                if len(fields) > 1:
                    self.datasets[fields[1]] = fields[2:]
        else:
            message("Database not found.")

    def query(self, query):
        message("Sending query", type="DEBUG")
        self.response = requests.get(self.url,
                                     query,
                                     proxies=self.proxies)

        message("Checking http response", type="DEBUG")

        if self.response.status_code != requests.codes.ok:
            msg = "HTTP response status code: {c}. {m}"
            msg = msg.format(c=str(self.response.status_code),
                             m=self.response.reason)
            message(msg, type="ERROR")

        msg = "([ \.\w]+ service you requested is currently unavailable[ \.\w]+)"
        hit = re.search(msg, self.response.text)

        if hit:
            msg = re.search(msg, self.response.text).group(1)
            message(msg.lstrip().rstrip(), type="WARNING")
            message("More information about this downtime "
                    "may be available on http://www.ensembl.info/",
                    type="ERROR")
