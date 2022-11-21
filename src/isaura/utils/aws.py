"""AWS utils."""

from decimal import Decimal
from typing import Any

import boto3


def get_dynamo_table(table_name: str) -> Any:
    """Create a boto3 dynamo client.

    Args:
        table_name (str): DynamoDB table name

    Returns:
        Any: DynamoDB table client
    """
    dynamo = boto3.resource("dynamodb", region_name="us-east-1")
    return dynamo.Table(table_name)


def decimal_json_decoder(obj):
    if isinstance(obj, Decimal):
        return float(obj)
    raise TypeError
