from typing import Any, Dict

from aws_lambda_powertools.logging import Logger
from aws_lambda_powertools.utilities.parser import parse
from aws_lambda_powertools.utilities.typing import LambdaContext
from pydantic import ValidationError

from isaura.blocs.precalc import read_precalc
from isaura.utils import get_dynamo_table


# NOTE: Relative import is important due to how lambdas work
# !WARN: Do not change to absolute imports
from .schema import ResponseSchema, RequestSchema


# Setup handler state
logger = Logger(service="isaura_precalc_get", level="INFO")
isaura_table = get_dynamo_table("isaura")


# API handler
def get(event: Dict[str, Any], context: LambdaContext) -> ResponseSchema:
    """Get handler."""
    try:
        req = parse(event, RequestSchema)
    except ValidationError as e:
        logger.exception("Validation Error")
        return ResponseSchema(e)

    query_params = req.queryStringParameters
    token = req.requestContext.authorizer.jwt

    try:
        precalc_list = read_precalc(isaura_table, token, query_params)
    except BaseException as e:
        # TODO: log error to a central service
        logger.exception("Query resolve error")
        return ResponseSchema(e)

    return ResponseSchema(precalc_list)
