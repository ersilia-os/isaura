from http import HTTPStatus
from typing import Any, Dict

from aws_lambda_powertools.logging import Logger
from aws_lambda_powertools.utilities.parser import parse
from aws_lambda_powertools.utilities.typing import LambdaContext
from pydantic import ValidationError

from isaura.blocs.precalc import read_precalc


# NOTE: Relative import is important due to how lambdas work
# !WARN: Do not change to absolute imports
from .schema import ResponseSchema, ResponseBodySchema, RequestSchema


# Setup handler state
logger = Logger(service="isaura_precalc_get", level="INFO")


# API handler
def get(event: Dict[str, Any], context: LambdaContext) -> ResponseSchema:
    """Get handler.

    Args:
        event (Dict[str, Any]): Lamda event
        context (LambdaContext): Lambda context

    Returns:
        ResponseSchema: API response
    """
    try:
        req = parse(event, RequestSchema)
    except ValidationError as e:
        logger.exception("Validation Error")
        return ResponseSchema(
            statusCode=HTTPStatus.BAD_REQUEST,
            headers={"Content-Type": "application/json"},
            body=ResponseBodySchema(
                msg="FAILED", items=[], last_eval_key=None, errors=e
            ),
        ).dict()

    query_params = req.queryStringParameters
    try:
        precalc_list, last_eval_key = read_precalc(query_params)
    except BaseException as e:
        # TODO: log error to a central service
        logger.exception("Query resolve error")
        return ResponseSchema(
            statusCode=HTTPStatus.BAD_REQUEST,
            headers={"Content-Type": "application/json"},
            body=ResponseBodySchema(
                msg="FAILED", items=[], last_eval_key=None, errors=e
            ),
        ).dict()

    return ResponseSchema(
        statusCode=HTTPStatus.OK,
        headers={"Content-Type": "application/json"},
        body=ResponseBodySchema(
            msg="OK", items=precalc_list, last_eval_key=last_eval_key, errors=[]
        ),
    ).dict()
