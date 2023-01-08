"""Isaura clients."""

import json
from pathlib import Path
from typing import List, Optional

import requests
from pydantic import parse_obj_as
from sqlalchemy import create_engine, delete, select
from sqlalchemy.orm import sessionmaker

from isaura.local.mixins import DBModel
from isaura.local.models import PrecalcDB
from isaura.routes.precalc.get.app.schema import (
    ResponseBodySchema,
    QueryParams,
    QueryType,
)
from isaura.routes.schemas.common import Precalc
from isaura.utils.aws import get_dynamo_table
from isaura.utils.dirs import get_workspace_path

default_workspace_path = get_workspace_path()


class IsauraAdminClient:
    def __init__(self: "IsauraAdminClient", table_name: str = "isaura_table") -> None:
        self.dynamo_table_client = get_dynamo_table(table_name)

    def insert(self: "IsauraAdminClient", precalc_list: List[Precalc]) -> None:
        """Add precals to dynamo table in batches.

        Args:
            precalc_list (List[Precalc]): List of Precalc objects
        """
        with self.dynamo_table_client.batch_writer() as batch:
            for precalc in precalc_list:
                batch.put_item(Item=precalc.dict(by_alias=True))

    def delete(self: "IsauraAdminClient", precalc_id_list: List[str]) -> None:
        """Remove precals from dynamo table in batches.

        Args:
            precalc_id_list (List[Precalc]): List of Precalc Ids to delete.
        """
        with self.dynamo_table_client.batch_writer() as batch:
            for precalc_id in precalc_id_list:
                batch.delete_item(
                    Key={"pk": precalc_id.split("#")[0], "sk": precalc_id.split("#")[1]}
                )


class IsauraRemoteClient:
    def __init__(self: "IsauraRemoteClient", api_url: str) -> None:
        self.api_url = api_url

    def get_all_precalcs(
        self: "IsauraRemoteClient", last_eval_key: Optional[str] = None
    ) -> ResponseBodySchema:
        """Get all precals.

        Args:
            last_eval_key (Optional[str], optional): Last eval key for pagination. Defaults to None.

        Returns:
            ResponseBodySchema: Api response.
        """
        params = QueryParams(
            query_type=QueryType.GET_ALL_PRECALC, last_eval_key=last_eval_key
        )

        resp = requests.get(url=self.api_url, params=params.dict())
        return ResponseBodySchema(
            msg=resp.json()["msg"],
            items=parse_obj_as(List[Precalc], resp.json()["items"]),
            last_eval_key=resp.json()["last_eval_key"],
            errors=resp.json()["errors"],
        )

    def get_precalc_by_id(
        self: "IsauraRemoteClient", precalc_id: str, only_keys: bool = False
    ) -> ResponseBodySchema:
        """Get precalc by id.

        Args:
            precalc_id (str): model_id#input_key.
            only_keys (bool): Only fetch keys if present.

        Returns:
            ResponseBodySchema: Api response.
        """
        params = QueryParams(
            query_type=QueryType.GET_PRECALC_BY_ID,
            last_eval_key=None,
            precalc_id=precalc_id,
            only_keys=only_keys,
        )

        resp = requests.get(url=self.api_url, params=params.dict())
        return ResponseBodySchema(
            msg=resp.json()["msg"],
            items=parse_obj_as(List[Precalc], resp.json()["items"]),
            last_eval_key=resp.json()["last_eval_key"],
            errors=resp.json()["errors"],
        )

    def get_precalcs_by_model_id(
        self: "IsauraRemoteClient", model_id: str
    ) -> ResponseBodySchema:
        """Get all precals.

        Args:
            model_id (str): model ID.

        Returns:
            ResponseBodySchema: Api response.
        """
        params = QueryParams(
            query_type=QueryType.GET_PRECALC_BY_MODEL_ID,
            last_eval_key=None,
            model_id=model_id,
        )

        resp = requests.get(url=self.api_url, params=params.dict())
        return ResponseBodySchema(
            msg=resp.json()["msg"],
            items=parse_obj_as(List[Precalc], resp.json()["items"]),
            last_eval_key=resp.json()["last_eval_key"],
            errors=resp.json()["errors"],
        )

    def get_precalcs_by_input_key(
        self: "IsauraRemoteClient", model_id_list: List[str], input_key: str
    ) -> ResponseBodySchema:
        """Get all precals.

        Args:
            model_id_list (str): List of model IDs.
            input_key (str): ID of the molecule.

        Returns:
            ResponseBodySchema: Api response.
        """
        items = []
        for model_id in model_id_list:
            resp = self.get_precalc_by_id(f"{model_id}#{input_key}")
            items.append(resp.items[0])

        return ResponseBodySchema(
            msg=resp.msg,
            items=items,
            last_eval_key=None,
            errors=resp.errors,
        )


class IsauraLocalClient:
    """Client to interact with local cache."""

    def __init__(
        self, db_path: Path = default_workspace_path.joinpath("isaura_local.db")
    ) -> None:
        self.db_engine = create_engine(f"sqlite:///{db_path}")
        self.DBSession = sessionmaker(self.db_engine)
        DBModel.metadata.create_all(self.db_engine)

    def reset(self: "IsauraLocalClient") -> None:
        """Reset the local cache database."""
        DBModel.metadata.drop_all(self.db_engine)
        DBModel.metadata.create_all(self.db_engine)

    def insert(self: "IsauraLocalClient", precalc_list: List[Precalc]) -> None:
        """Add precals to local sqlitedb in batches.

        Args:
            precalc_list (List[Precalc]): List of Precalc objects
        """
        with self.DBSession() as session:
            precalc_db_list = [
                PrecalcDB(precalc_id=precalc.precalc_id, **json.loads(precalc.json()))
                for precalc in precalc_list
            ]
            PrecalcDB.bulk_create(precalc_db_list, session)

    def delete(self: "IsauraLocalClient", precalc_id_list: List[str]) -> None:
        """Remove precals from local sqlitedb in batches.

        Args:
            precalc_id_list (List[Precalc]): List of Precalc Ids to delete.
        """
        stmt = delete(PrecalcDB).where(PrecalcDB.precalc_id.in_(precalc_id_list))
        print(stmt)
        with self.DBSession() as session:
            session.execute(stmt)
            session.commit()

    def get_all_precalcs(
        self: "IsauraLocalClient", page: int = 0, limit: int = 100
    ) -> List[Precalc]:
        """Get all precals.

        Args:
            page (int): page key for pagination. Defaults to 0.
            limit (int): Items to return per page.

        Returns:
            List[Precalc]: List of fetched precalcs.
        """
        stmt = select(PrecalcDB).limit(limit).offset(page * limit)
        with self.DBSession() as session:
            precalc_db_list: List[PrecalcDB] = session.execute(stmt).scalars().all()
            return [Precalc.from_orm(obj) for obj in precalc_db_list]

    def get_precalc_by_id(self: "IsauraLocalClient", precalc_id: str) -> List[Precalc]:
        """Get precalc by id.

        Args:
            precalc_id (str): model_id#input_key.

        Returns:
            List[Precalc]: List of fetched precalcs.
        """
        stmt = select(PrecalcDB).where(PrecalcDB.precalc_id == precalc_id)

        with self.DBSession() as session:
            precalc_db_list: List[PrecalcDB] = session.execute(stmt).scalars().all()
            return [Precalc.from_orm(obj) for obj in precalc_db_list]

    def get_precalc_by_model_id(
        self: "IsauraLocalClient", model_id: str, page: int = 0, limit: int = 100
    ) -> List[Precalc]:
        """Get precalc by id.

        Args:
            model_id (str): model_id.
            page (int): page key for pagination. Defaults to 0.
            limit (int): Items to return per page.

        Returns:
            List[Precalc]: List of fetched precalcs.
        """
        stmt = (
            select(PrecalcDB)
            .where(PrecalcDB.model_id == model_id)
            .limit(limit)
            .offset(page * limit)
        )

        with self.DBSession() as session:
            precalc_db_list: List[PrecalcDB] = session.execute(stmt).scalars().all()
            return [Precalc.from_orm(obj) for obj in precalc_db_list]

    def get_precalc_by_input_key(
        self: "IsauraLocalClient", input_key: str, page: int = 0, limit: int = 100
    ) -> List[Precalc]:
        """Get precalc by id.

        Args:
            input_key (str): input_key.
            page (int): page key for pagination. Defaults to 0.
            limit (int): Items to return per page.

        Returns:
            List[Precalc]: List of fetched precalcs.
        """
        stmt = (
            select(PrecalcDB)
            .where(PrecalcDB.input_key == input_key)
            .limit(limit)
            .offset(page * limit)
        )

        with self.DBSession() as session:
            precalc_db_list: List[PrecalcDB] = session.execute(stmt).scalars().all()
            return [Precalc.from_orm(obj) for obj in precalc_db_list]
