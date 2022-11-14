"""Adjustment GET API schemas."""

from enum import Enum
from typing import Any, List, Optional

from pydantic import BaseModel

from isaura.routes.schemas.common import RequestContext


class Precalc(BaseModel):
    """Precalulated predictions"""

    pass


class FilterParams(BaseModel):
    """Filter parameters."""

    precalc_id: str


class OrderByEnum(Enum):
    """OrderBy Enum."""

    CREATION_TIMESTAMP = "CREATION_TIMESTAMP"


class QueryType(Enum):
    """Query type enum."""

    GET_ALL_PRECALC = "get_all_precalc"
    GET_PRECALC_BY_ID = "get_precalc_by_id"


class QueryParams(BaseModel):
    """Query parameters."""

    query_type: QueryType
    last_eval_key: Optional[int]
    precalc_id: Optional[str]


class RequestSchema(BaseModel):
    """Request schema."""

    queryStringParameters: QueryParams
    requestContext: RequestContext


class ResponseBodySchema(BaseModel):
    """Reponse body schema."""

    msg: str
    items: List[Precalc]
    errors: Any


class ResponseSchema(BaseModel):
    """Response schema."""

    statusCode: int
    headers: Any
    body: ResponseBodySchema
