"""AWS Precalculation storage infrastructure stack."""

from pathlib import Path
from os import path, makedirs
import shutil

from aws_cdk import (
    aws_dynamodb as dynamodb,
    aws_lambda as λ,
    aws_logs as logs,
    Stack,
    CfnOutput,
)

from aws_cdk.aws_s3 import Bucket, HttpMethods as S3HttpMethod


class IsauraMainStack(Stack):
    """Isaura main stack."""

    def create_lambda_layer(self) -> None:
        """Create packages zip."""
        # Create di
        self.project_root = Path.cwd()
        makedirs(
            self.project_root.joinpath(".venv_prod/layer/packages/python"),
            exist_ok=True,
        )

        # Copy site packages
        shutil.copytree(
            src=self.project_root.joinpath(
                ".venv_prod/lib/python3.8/site-packages"
            ).absolute(),
            dst=self.project_root.joinpath(
                ".venv_prod/layer/packages/python"
            ).absolute(),
            dirs_exist_ok=True,
        )

        # Create zip file
        shutil.make_archive(
            base_name=self.project_root.joinpath(".venv_prod/layer/python"),
            format="zip",
            root_dir=self.project_root.joinpath(".venv_prod/layer/packages"),
        )

    def __init__(self: "IsauraMainStack", scope: Stack, id: str, **kwargs) -> None:
        """Initialize Bob's Engine main stack.

        Args:
            scope (core.Stack): Scope of the stack
            id (str): Stack ID
        """
        super().__init__(scope, id, **kwargs)

        # S3 bucket
        isaura_bucket = Bucket(
            scope=self, id="isaura_bucket", bucket_name="isaura-bucket"
        )

        isaura_bucket.add_cors_rule(
            allowed_methods=[
                S3HttpMethod.GET,
                S3HttpMethod.POST,
                S3HttpMethod.PUT,
                S3HttpMethod.DELETE,
                S3HttpMethod.HEAD,
            ],
            allowed_origins=[
                "http://localhost",
                "http://localhost:8888",
                "https://amazonaws.com",
            ],
            allowed_headers=["*"],
        )

        # DynamoDB
        isaura_table = dynamodb.Table(
            scope=self,
            id="isaura_table",
            table_name="isaura_table",
            billing_mode=dynamodb.BillingMode.PAY_PER_REQUEST,
            partition_key=dynamodb.Attribute(
                name="pk", type=dynamodb.AttributeType.STRING
            ),
            sort_key=dynamodb.Attribute(name="sk", type=dynamodb.AttributeType.STRING),
        )

        # Add GSI to the table
        # isaura_table.add_global_secondary_index()

        # Create lambda layer
        self.create_lambda_layer()

        # Shared lambda layer
        isaura_lambda_layer = λ.LayerVersion(
            scope=self,
            id="isaura_lambda_layer",
            code=λ.Code.from_asset(
                path.abspath(
                    self.project_root.joinpath(".venv_prod/layer/python.zip").absolute()
                )
            ),
            compatible_runtimes=[λ.Runtime.PYTHON_3_8],
            description="Lambda layer for isaura functions",
        )
        isaura_lambda_layer.add_permission("account-grant", account_id=self.account)

        # Lambda funtions

        # Status
        status_heartbeat_λ = λ.Function(
            scope=self,
            id="isaura_lambda_status_heartbeat",
            function_name="isaura_lambda_status_heartbeat",
            runtime=λ.Runtime.PYTHON_3_8,
            handler="app.handler.get",
            layers=[isaura_lambda_layer],
            code=λ.Code.from_asset(
                path.abspath(
                    Path(__file__).parents[1].joinpath("routes/status/get").absolute()
                )
            ),
            log_retention=logs.RetentionDays.THREE_DAYS,
        )
        status_heartbeat_λ_url = status_heartbeat_λ.add_function_url(
            auth_type=λ.FunctionUrlAuthType.NONE
        )
        isaura_table.grant_read_data(status_heartbeat_λ)

        precalc_get_λ = λ.Function(
            scope=self,
            id="isaura_lambda_precalc_get",
            function_name="isaura_lambda_precalc_get",
            runtime=λ.Runtime.PYTHON_3_8,
            handler="app.handler.get",
            layers=[isaura_lambda_layer],
            code=λ.Code.from_asset(
                path.abspath(
                    Path(__file__).parents[1].joinpath("routes/precalc/get").absolute()
                )
            ),
            log_retention=logs.RetentionDays.THREE_DAYS,
        )
        precalc_get_λ_url = precalc_get_λ.add_function_url(
            auth_type=λ.FunctionUrlAuthType.NONE
        )
        isaura_table.grant_read_data(precalc_get_λ)

        CfnOutput(self, "region", value=self.region)
        CfnOutput(self, "statusLambdaUrl", value=status_heartbeat_λ_url.url)
        CfnOutput(self, "precalcGetLambdaUrl", value=precalc_get_λ_url.url)
