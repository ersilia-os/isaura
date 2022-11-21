"""Common schemas."""

from decimal import Decimal
import json
from typing import Dict

from pydantic import BaseModel, Field, validator

from isaura.utils.aws import decimal_json_decoder


class RequestContext(BaseModel):
    """Common request context model."""

    pass


class Precalc(BaseModel):
    """Schema for Preclac object."""

    model_id: str = Field(alias="pk")
    input_key: str = Field(alias="sk")
    value: Dict

    @validator("value")
    def convert_float_to_decimals(cls, v):
        return json.loads(
            json.dumps(v, default=decimal_json_decoder), parse_float=Decimal
        )

    @property
    def precalc_id(self: "Precalc") -> str:
        """Get precalc id str.

        Returns:
            str: precalc_id
        """
        return f"{self.model_id}#{self.input_key}"

    class Config:
        ignore_extra = False
        allow_population_by_field_name = True
