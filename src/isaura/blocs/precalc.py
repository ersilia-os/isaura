"""Precalc BLoCs."""

from typing import List, Tuple

from boto3.dynamodb.conditions import Key
from pydantic import parse_obj_as, ValidationError

from isaura.routes.precalc.get.app.schema import QueryParams, QueryType
from isaura.routes.schemas.common import Precalc
from isaura.utils.aws import get_dynamo_table


def read_precalc(query_params: QueryParams) -> Tuple[List[Precalc], str]:
    """Read precals from the dynamo table.

    Args:
        query_params (QueryParams): Query parameters for dynamo search

    Raises:
        e : Query params validation error based on context of search

    Returns:
        Tuple[List[Precalc], str] : A tuple of precalc list and optional last eval key
    """
    dynamo_table_client = get_dynamo_table("isaura_table")
    last_eval_key = None
    res = None
    query_type = QueryType(query_params.query_type)

    if query_type == QueryType.GET_ALL_PRECALC:
        res = dynamo_table_client.scan()

    elif query_type == QueryType.GET_PRECALC_BY_ID:
        if query_params.precalc_id is None:
            e = ValidationError(
                errors=[
                    {
                        "loc": "precalc_id",
                        "msg": "field required",
                        "type": "value_error.missing",
                    }
                ],
                model=QueryParams,
            )
            raise e
        res = dynamo_table_client.query(
            KeyConditionExpression=Key("pk").eq(query_params.precalc_id.split("#")[0])
            & Key("sk").eq(query_params.precalc_id.split("#")[1])
        )
    elif query_type == QueryType.GET_PRECALC_BY_MODEL_ID:
        if query_params.model_id is None:
            e = ValidationError(
                errors=[
                    {
                        "loc": "model_id",
                        "msg": "field required",
                        "type": "value_error.missing",
                    }
                ],
                model=QueryParams,
            )
            raise e
        res = dynamo_table_client.query(
            KeyConditionExpression=Key("pk").eq(query_params.model_id)
        )
    elif query_type == QueryType.GET_PRECALC_BY_INPUT_KEY:
        if query_params.model_id is None:
            e = ValidationError(
                errors=[
                    {
                        "loc": "model_id",
                        "msg": "field required",
                        "type": "value_error.missing",
                    }
                ],
                model=QueryParams,
            )
            raise e

        if query_params.input_key is None:
            e = ValidationError(
                errors=[
                    {
                        "loc": "input_key",
                        "msg": "field required",
                        "type": "value_error.missing",
                    }
                ],
                model=QueryParams,
            )
            raise e

        res = dynamo_table_client.query(
            KeyConditionExpression=Key("pk").eq(query_params.model_id)
            & Key("sk").eq(query_params.input_key)
        )
    if res is not None:
        res_items = parse_obj_as(List[Precalc], res["Items"])
        if "LastEvaluatedKey" in res:
            last_eval_key = res["LastEvaluatedKey"]
    else:
        res_items = []
    return res_items, last_eval_key


def write_precalc(precalc_list: List[Precalc]) -> List[Precalc]:
    pass


def update_precalc(precalc_list: List[Precalc]) -> List[Precalc]:
    pass


def delete_precalc(precalc_id_list: List[str]) -> None:
    pass
