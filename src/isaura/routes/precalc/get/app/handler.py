import json
from typing import Any, Dict, List, Optional

from aws_lambda_powertools.logging import Logger
from aws_lambda_powertools.utilities.parser import parse
from aws_lambda_powertools.utilities.typing import LambdaContext
import boto3
from boto3.dynamodb.conditions import Key
import botocore
from pydantic import parse_obj_as, ValidationError
from rivets.schemas.adjustment import Adjustment
from rivets.schemas.common import Config, Token
from rivets.utils import decimal_json_decoder


# NOTE: Relative import is important due to how lambdas work
# !WARN: Do not change to absolute imports
from .schema import OrderByEnum, QueryType, ResponseSchema, RequestSchema, QueryParams

logger = Logger(service="pahal_sre_bob_adjustment_get", level="INFO")


# API handler
def get(event: Dict[str, Any], context: LambdaContext) -> ResponseSchema:
    """Get handler."""
    try:
        req = parse(event, RequestSchema)
    except ValidationError as e:
        # TODO: log error to a central service
        logger.exception("Validation Error")
        return validation_error_response(e)

    query_params = req.queryStringParameters
    token = req.requestContext.authorizer.jwt

    try:
        precalc = read_precalc(isaura_table, token, query_params)
    except BaseException as e:
        # TODO: log error to a central service
        logger.exception("Query resolve error")
        return validation_error_response(e)

    return success_response("OK", [precalc])