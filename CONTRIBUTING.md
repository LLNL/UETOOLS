# Contributing to UETOOLS

Thank you for your interest in contributing to **UETOOLS**!
We welcome contributions from both core developers and the wider open-source community.

This document outlines the development workflow, branching strategy, and expectations for contributions.

---

## Code of Conduct

All contributors are expected to interact respectfully and professionally.
Harassment, discrimination, or unconstructive behavior will not be tolerated.

---

## Repository Overview

UETOOLS follows a **GitFlow-like** branching model designed to:

- Keep the `main` branch mostly stable
- Support active development on `develop`
- Enable clean releases and urgent hotfixes
- Encourage reviewable, high-quality contributions

---

## Branching Model

### Long-lived branches

#### `main`
- Represents the latest released or release-ready state
- Mostly stable
- Protected: **no direct commits allowed**
- Releases are tagged from this branch

#### `develop`
- Integration branch for ongoing development
- All features, bug fixes, and hotfixes eventually merge here
- Protected: **no direct commits allowed**

---

### Short-lived branches

| Branch type | Purpose | Base branch |
|------------|--------|-------------|
| `feature/*` | New features or enhancements | `develop` |
| `bugfix/*`  | Non-urgent bug fixes | `develop` |
| `hotfix/*`  | Urgent fixes for released versions | `main` |
| `release/*` | (Optional) Release stabilization | `develop` |

---

## Development Workflow

All changes must be submitted via Pull Requests.
Direct commits to `main` or `develop` are not allowed.

---

## Releases & Versioning

- Releases are tagged directly from `main`
- Semantic versioning is used: `vX.Y.Z`

---

## Issues

We use **GitHub Issues** to track bugs, feature requests, and discussions.

---

## Notes for External Contributors

External contributions are encouraged.

- Fork the repository
- Create branches following conventions
- Open PRs against `develop`
