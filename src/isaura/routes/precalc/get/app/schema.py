"""Adjustment GET API schemas."""

from enum import Enum
from typing import Any, List, Optional

from pydantic import BaseModel, validator

from isaura.routes.schemas.common import Precalc


class QueryType(Enum):
    """Query type enum."""

    GET_ALL_PRECALC = "get_all_precalc"
    GET_PRECALC_BY_ID = "get_precalc_by_id"
    GET_PRECALC_BY_MODEL_ID = "get_precalc_by_model_id"
    GET_PRECALC_BY_INPUT_KEY = "get_precalc_by_input_key"


class QueryParams(BaseModel):
    """Query parameters."""

    query_type: QueryType
    last_eval_key: Optional[int]
    precalc_id: Optional[str]
    model_id: Optional[str]
    input_key: Optional[str]

    @validator("query_type")
    def convert_query_type_to_str(cls, v):
        return v.value


class RequestSchema(BaseModel):
    """Request schema."""

    queryStringParameters: QueryParams


class ResponseBodySchema(BaseModel):
    """Reponse body schema."""

    msg: str
    items: List[Precalc]
    last_eval_key: Optional[str]
    errors: Any


class ResponseSchema(BaseModel):
    """Response schema."""

    statusCode: int
    headers: Any
    body: ResponseBodySchema
