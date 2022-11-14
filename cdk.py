#!/usr/bin/env python3

from aws_cdk import App

from isaura.iac.aws_precalc_stack import IsauraMainStack


app = App()

IsauraMainStack(app, "isaura")

app.synth()
