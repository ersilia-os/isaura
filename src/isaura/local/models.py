"""Sqlalchemy models."""

from sqlalchemy import Column, String, JSON

from isaura.local.mixins import DBModel


class PrecalcDB(DBModel):
    """Precalc ORM class."""

    __tablename__ = "precalc"
    precalc_id = Column(String, nullable=False, index=True, primary_key=True)
    model_id = Column(String, nullable=False, index=True)
    input_key = Column(String, nullable=False, index=True)
    value = Column(JSON, nullable=False)
