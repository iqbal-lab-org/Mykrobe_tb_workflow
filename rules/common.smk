from typing import List
import re


class InvalidBarcode(Exception):
    __module__ = Exception.__module__


def barcode_parser(barcodes_string: str) -> List[str]:
    """Parses the barcodes string and ensures they follow correct format"""
    msg = "Barcode must be of the form BC01. That is, BC followed by 2 digits."
    regex = r"\bBC\d{2}\b"
    barcodes = barcodes_string.split()
    for barcode in barcodes:
        if not (len(barcode) == 4 and re.match(regex, barcode)):
            raise InvalidBarcode(barcode + "\n" + msg)
    return barcodes


def infer_barcode_dir(wildcards) -> str:
    barcode_num = wildcards.sample[-2:]
    return "barcode{}".format(barcode_num)
